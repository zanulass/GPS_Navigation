module compute_receiver_position
  use WGS84_constants
  use exec_conditions
  use time_util
  use print_util
  use navigation_message
  use correction_parameters
  use compute_solution
  implicit none
contains
  subroutine main_calc_position(wt, pseudo_range, sol)
    implicit none
    ! 引数詳細
    TYPE(wtime), INTENT(INOUT)       :: wt
    DOUBLE PRECISION, INTENT(INOUT)  :: pseudo_range(MAX_PRN)
    DOUBLE PRECISION, INTENT(INOUT)  :: sol(MAX_UNKNOWNS)

    ! 使用局所領域
    INTEGER            :: return_code
    INTEGER            :: prn
    INTEGER            :: loop ! ループカウンタ
    DOUBLE PRECISION   :: clock_err  ! クロック誤差(秒に変換)


    ! エフェメリスをセット
    do prn=1, MAX_PRN
      if (pseudo_range(prn) > 0.d0) then
        ! 観測量が設定された衛星に対して，エフェメリス情報をセット
        call set_ephemeris(prn, wt, -1, return_code)
        ! セットしたエフェメリス情報を実行結果リストに出力
        call print_ephemeris_info(prn)

        ! 有効期限内のエフェメリス情報がなければ無効な衛星として擬似距離に0を設定
        if (return_code == 9) then
          pseudo_range(prn) = 0.d0
        ! healthフラグをチェックし，無効な衛星は擬似距離に0を設定
        else if (current_ephem(prn)%health /= 0.d0) then
          pseudo_range(prn) = 0.d0
        end if
      end if
    end do

    ! 解を初期化
    sol(:) = 0.d0
    ! sol(1) = -3947762.7496d0
    ! sol(2) = 3364399.8789d0
    ! sol(3) = 3699428.5111d0
    ! sol(4) = 0.d0

    ! 測位計算開始
    do loop=1, MAX_LOOP
      ! 受信機位置計算
      call calc_position(wt, pseudo_range, sol, loop, return_code)
      if (return_code /= 0) exit

      ! 受信機クロック誤差
      clock_err = sol(4) / C

      ! 推定パラメータを保存
      sol_for_print(1:3, loop) = sol(1:3)
      sol_for_print(4, loop) = clock_err

      ! 途中経過を出力
      write(6, '(A,I0, 5X, A, f12.3,5X, A ,f12.3,5X, A, f12.3, 5X, A, D12.3)') &
      "LOOP = ", loop, "x = ", sol(1), "y = ", sol(2), "z = ", sol(3), "s = ", clock_err

    end do

  end subroutine main_calc_position

  subroutine calc_position(wt, pseudo_range, sol, loop, return_code)
    implicit none
    ! 引数詳細
    TYPE(wtime), INTENT(INOUT)       :: wt
    DOUBLE PRECISION, INTENT(INOUT)  :: pseudo_range(MAX_PRN)
    DOUBLE PRECISION, INTENT(INOUT)  :: sol(MAX_UNKNOWNS)
    INTEGER, INTENT(IN)              :: loop  ! csv書き出し用にループ番号を取得
    INTEGER, INTENT(INOUT)           :: return_code


    ! 使用局所領域
    DOUBLE PRECISION  :: delta_pseudo_range(MAX_PRN)  ! 擬似距離の残差
    DOUBLE PRECISION  :: iono_correction  ! 電離層遅延補正値
    DOUBLE PRECISION  :: delta_pseudo_range_iono(MAX_PRN)  ! 擬似距離の残差(電離層遅延成分)
    DOUBLE PRECISION  :: tropo_correction  ! 対流圏遅延補正値
    DOUBLE PRECISION  :: delta_pseudo_range_tropo(MAX_PRN)  ! 擬似距離の残差(対流圏遅延成分)
    INTEGER           :: prn
    DOUBLE PRECISION  :: range_correct_clock  ! クロック補正値計算用のrange
    DOUBLE PRECISION  :: range_calc_sat_pos  ! 衛星位置計算用のrange
    DOUBLE PRECISION  :: euclidian_distance  ! 観測行列に使用するrange
    DOUBLE PRECISION  :: sat_pos(3)  ! ECEF座標系における衛星位置ベクトル
    DOUBLE PRECISION  :: receiver_pos(3)  ! ECEF座標系における受信機位置ベクトル
    DOUBLE PRECISION  :: sat_clock  ! クロック補正値
    DOUBLE PRECISION  :: obs_mat(MAX_SATS, MAX_UNKNOWNS)  ! 観測行列(観測衛星数の上限 ×　未知数の上限)
    DOUBLE PRECISION  :: delta_range(MAX_SATS)  ! rangeの修正量
    DOUBLE PRECISION  :: delta_x(MAX_UNKNOWNS)  ! 解の更新量
    INTEGER           :: rtn = 0 ! 最小二乗法サブルーチンのリターンコード
    INTEGER           :: i
    INTEGER           :: u ! デバッグ用

    ! リターンコードを初期化
    return_code = 0

    ! 使用済み衛星数を初期化
    num_used_PRN = 0
    ! 測位計算に使用する衛星のPRNを格納するリストを初期化
    used_PRN_list(:) = 0

    ! 受信時刻(受信機クロック分を補正する)
    wt%sec = wt%sec - sol(4) / C

    ! 暫定の受信機位置
    receiver_pos(1:3) = sol(1:3)

    !  観測行列を初期化
    obs_mat(:, :) = 0.d0

    do prn=1, MAX_PRN
      ! 残差を初期化
      delta_pseudo_range(prn) = 0.d0
      delta_pseudo_range_iono(prn) = 0.d0
      delta_pseudo_range_tropo(prn) = 0.d0
      ! 擬似距離が有効な衛星のみを使用する(擬似距離が0のPRNは飛ばす)
      if (pseudo_range(prn) == 0.d0) cycle

      ! *** 有効な衛星のPRNを記録 ***
      num_used_PRN = num_used_PRN + 1
      used_PRN_list(num_used_PRN) = prn

      ! *** 衛星クロック補正値を計算 ***
      ! 衛星クロック補正値計算用のrangeを設定
      range_correct_clock = pseudo_range(prn) - sol(4)
      sat_clock = 0.d0
      call correct_sat_clock(prn, wt, range_correct_clock, sat_clock)

      ! *** 衛星位置を計算 ***
      ! 衛星位置計算用のrangeを設定
      range_calc_sat_pos = range_correct_clock + (sat_clock * C)
      sat_pos(:) = 0.d0
      call calc_sat_position(prn, wt, range_calc_sat_pos, sat_pos)

      ! *** 観測行列(observastion matrix)の作成 ***
      euclidian_distance = sqrt( sum( (sat_pos(1:3) - receiver_pos(1:3)) ** 2.d0 ) )
      obs_mat(num_used_PRN, 1:3) = ( sol(1:3) - sat_pos(1:3) ) / euclidian_distance
      obs_mat(num_used_PRN, 4) = 1.d0

      ! *** 擬似距離の修正量を計算 ***
      delta_range(num_used_PRN) = pseudo_range(prn) + (sat_clock * C) &
                                    - (euclidian_distance + sol(4) )

      ! *** 電離層補正補正 ***
      iono_correction = 0.d0
      if (iono_flag .eqv. .true.) then
        call calc_iono_correction(sat_pos, receiver_pos, wt, iono_correction)
        delta_pseudo_range_iono(prn) = delta_pseudo_range_iono(prn) + iono_correction

        ! 擬似距離の修正量に電離層遅延補正量を加える
        delta_range(num_used_PRN) = delta_range(num_used_PRN) + delta_pseudo_range_iono(prn)
      end if

      ! *** 対流圏遅延補正 ***
      tropo_correction = 0.d0
      if (tropo_flag .eqv. .true.) then
        call calc_tropo_correction(sat_pos, receiver_pos, tropo_correction)
        delta_pseudo_range_tropo(prn) = delta_pseudo_range_tropo(prn) + tropo_correction

        ! 擬似距離の修正量に対流圏遅延補正量を加える
        delta_range(num_used_PRN) = delta_range(num_used_PRN) + delta_pseudo_range_tropo(prn)
      end if

      ! *** 出力csv用のデータを保存 ***
      ! 擬似距離の残差を保存
      oc_for_print(prn, loop) = delta_range(num_used_PRN)

      ! 衛生位置を保存
      sat_pos_for_print(1:3, prn, loop) = sat_pos

      ! 補正値パラメータを保存
      sat_clock_for_print(prn, loop) = sat_clock
      iono_correction_for_print(prn, loop) = iono_correction
      tropo_correction_for_print(prn, loop) = tropo_correction


    end do

    ! 未知数(x, y, z, s)に対して有効な衛星数が足りなければエラー
    if (num_used_PRN < MAX_UNKNOWNS) then
      return_code = 9
      goto 8000
    end if

    ! 方程式を解く
    call least_squares(obs_mat, delta_range, delta_x, num_used_PRN, 4, rtn)
    if (rtn /= 0) then
      return_code = 9
      goto 9000
    end if

    ! 解の更新
    sol(1:4) = sol(1:4) + delta_x(1:4)






    return

    ! *** エラー処理 ***
    8000 continue
    write(*, *) "使用できる衛星数が未知数の数より少ない"
    return

    9000 continue
    write(*, *) "最小二乗法の計算に失敗"
    return

  end subroutine calc_position

end module compute_receiver_position