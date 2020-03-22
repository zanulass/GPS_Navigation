program main
  use mod_variable
  use time_util
  use navigation_message
  use satellite_position
  use print_list
  use compute_solution
  implicit none

!   1. 機能概要
!     GPS位置計算メインルーチン
!     RINEX形式の航法データと観測データを用いて、受信機の位置を計算する

!   2. 注意事項
!     特になし
!
!   3. 使用モジュール
!     main(本ルーチン)
!       +
!       +----mod_variable (共通の定数と変数を宣言)
!       +
!       +----time_util (衛星クロックの補正)
!       +
!       +----navigation_message (Navigation Message Fileの読み込み)
!       +
!       +----satellite_position (衛星の位置計算)
!       +
!       +----plint_list (エフェメリス情報と衛星位置を実行結果リストファイルに書き出し)
!       +
!       +----compute_solutionモジュール (最小二乗法により方程式を解く)
! 　
! 　4 引数詳細
!     なし

!   5. 使用局所領域
!   +------------------------------------------------------------------------------------------------------------
!   ! TYPE(LENGTH)            !  name(size)
!   +------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER      :: SATS = 5                    ! 使用する衛星数
      INTEGER, PARAMETER      :: MAX_LOOP = 8                ! 解を求める際に用いる最大ループ回数
      INTEGER                 :: PRN_list(SATS)              ! 使用する衛星のリスト
      DOUBLE PRECISION        :: SATS_RANGE(SATS)            ! 各衛星の位置とrangeの配列
      DOUBLE PRECISION        :: r                           ! range
      DOUBLE PRECISION        :: G(SATS, MAX_UNKNOWNS)   ! 観測行列(観測衛星数の上限 ×　未知数の上限)
      DOUBLE PRECISION        :: dr(SATS)                    ! rangeの修正量
      DOUBLE PRECISION        :: dx(SATS)                    ! 解の更新量
      DOUBLE PRECISION        :: sol(MAX_UNKNOWNS)           ! 方程式の解:受信機位置x,y,z座標, 受信機クロック誤差 s
      INTEGER                 :: i, n, loop, u               ! ループ用カウンタ
      CHARACTER(256)          :: nav_msg_file                ! RINEX NAVIGATION MESSAGE FILEのパス
      CHARACTER(256)          :: list_file                   ! 実行結果リストのパス
      TYPE(wtime)             :: wt                          ! 時刻
      TYPE(ephemeris_info)    :: current_ephem               ! 作業用のエフェメリス情報
      DOUBLE PRECISION        :: sat_clock                   ! 衛星のクロック補正量

      DOUBLE PRECISION        :: x, y, z, s                    ! 解の確認用出力
!   +-----------------------------------------------------------------------------------------------------------------


  ! SAT1_POS(:) = (/ -13897607.6294d0, -10930188.6233d0, 19676689.6804d0 /)
  ! SAT2_POS(:) = (/ -17800899.1998d0, 15689920.8120d0, 11943543.3888d0 /)
  ! SAT3_POS(:) = (/ -1510958.2282d0, 26280096.7818d0, -3117646.1949d0 /)
  ! SAT4_POS(:) = (/ -12210758.3517d0, 20413597.0201d0, -11649499.5474d0 /)
  ! SAT5_POS(:) = (/ -170032.6981d0, 17261822.6784d0, 20555984.4061d0 /)
  ! SATS_POSITION(:, :) = 0.0d0 ! 衛星の位置を初期化
  ! do i = 1, 3 ! 各衛星の位置を配列にセット
  !   SATS_POSITION(i, 1) = SAT1_POS(i)
  !   SATS_POSITION(i, 2) = SAT2_POS(i)
  !   SATS_POSITION(i, 3) = SAT3_POS(i)
  !   SATS_POSITION(i, 4) = SAT4_POS(i)
  !   SATS_POSITION(i, 5) = SAT5_POS(i)
  ! end do

  ! 使用するPRN
  PRN_list(:) = (/ 5,14,16,22,25 /)
  ! 観測データを配列にセット
  ! sats_range(:) = (/ 23634878.5219d0, 20292688.3557d0, 24032055.0372d0, 24383229.3740d0, 22170992.8178d0/)

  sats_range(:) = (/ &
    23545777.534d0, & ! PRN 05
    20299789.570d0, & ! PRN 14
    24027782.537d0, & ! PRN 16
    24367716.061d0, & ! PRN 22
    22169926.127d0  & ! PRN 25
   /)

  ! Navigatione Message Fileのパス
  nav_msg_file = "../data/mtka3180.05n"
  list_file = "../tmp/list"



  ! 時刻を指定
  wt%week = 1349     ! 05/11/13〜19の週
  wt%sec = 86400.d0  ! 月曜日の00:00:00



  open (20,  file=list_file, action='write', status='replace') ! 実行結果リストオープン

  call read_nav_msg(nav_msg_file) ! Navigation Message File読み込み
  call print_nav_file_header() ! Navigation Message ヘッダ部をリストに書き出し

  sol(:) = 0.d0  ! 解を初期化
  G(:, :) = 0.d0 ! 観測行列を初期化
  ! 解を求めるループ
  do loop=1, MAX_LOOP
    n = SATS  ! 衛星の数をセットする
    do i=1, n ! 衛星の数だけループする
      call set_ephemeris(PRN_list(i), wt, -1, current_ephem) ! 作業中の衛星のエフェメリスをセット

      sat_clock = 0.d0 ! 衛星のクロック補正量を初期化
      call correct_sat_clock(wt, current_ephem, sat_clock) ! Navigation Message ヘッダ部からクロック補正量sat_clockを計算

      current_ephem%pos_xyz(:) = 0.d0 ! 衛星位置を初期化
      call calc_satpos(wt, current_ephem) ! 衛星の位置を計算　ECEFにおける衛星座標がpos_xyzにセットされる

      if (loop == 8) then
        call print_ephemeris_info(current_ephem) ! 最後に一度、各衛星のエフェメリスをリストに書き出し
      end if

      r = sqrt( sum( (current_ephem%pos_xyz(1:3) - sol(1:3) ) ** 2.d0 ) ) ! 疑似距離rの計算
      G(i,1:3) = ( sol(1:3) - current_ephem%pos_xyz(1:3) ) / r ! 観測行列Gを作成
      G(i, 4) = 1.d0

      dr(i) = SATS_RANGE(i) + sat_clock*C - (r + sol(4)) !擬似距離の修正量drを計算
    end do
    ! 観測行列のデバッグライト
      ! write(*, *) '========================================='
      ! do u = 1, n
      !   write(*, *) G(u, 1:4)
      ! end do

    call least_squares(G, dr, dx, n, 4); ! 最小二乗法により方程式を解く

    do i=1, 4
      sol(i) = sol(i) + dx(i) ! 解を更新
    end do

    ! 途中経過を出力
    write(6, '("LOOP ",I0, 5X,"x = ",f12.3,5X,"y = ",f12.3,5X,"z = ", f12.3, 5X,"s = ", D12.3)') &
      loop, sol(1), sol(2), sol(3), sol(4)
  end do
  close(20) ! 実行結果リストクローズaaaaaa

  ! ! 正しい解
  x = -3947846.647d0
  y = 3364338.022d0
  z = 3699406.626d0
  s = -5.3233d-008
  write(6, *) "*****************************************************"
  write(6, '("x = ",f12.3,5X,"y = ",f12.3,5X,"z = ",f12.3, 5X,"s = ",d12.4)') x, y, z, s

end program main
