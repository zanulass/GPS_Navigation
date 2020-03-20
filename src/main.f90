program main
  use compute_solution
  implicit none

!   1. 機能概要
!     GPS位置計算メインルーチン
!     RINEX形式の航法データと観測データを用いて、受信機の位置を計算する

!   2. 注意事項
!     特になし

!   3. 外部効果
!     特になし

!   4. 使用モジュール
!     main(本ルーチン)
!       +
!       +----compute_solutionモジュール (最小二乗法により方程式を解く)
! 　
! 　5 引数詳細
!     なし

!   6. 使用局所領域
!   +------------------------------------------------------------------------------------------------------------------
!   ! TYPE(LENGTH)            !  name(size)
!   +-----------------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER      :: SATS = 5             ! 使用する衛星数
      INTEGER, PARAMETER      :: MAX_LOOP = 8         ! 解を求める際に用いる最大ループ回数
      DOUBLE PRECISION        :: SAT1_POS(3), SAT2_POS(3), SAT3_POS(3), SAT4_POS(3), SAT5_POS(3) ! 各衛星の位置
      DOUBLE PRECISION        :: SATS_POSITION(3, SATS), SATS_RANGE(SATS) ! 各衛生の位置とrangeの配列
      DOUBLE PRECISION        :: r                    ! range
      DOUBLE PRECISION        :: xyz(3)               ! 衛星位置xyzの配列
      DOUBLE PRECISION        :: G(SATS, 3)           ! 観測行列
      DOUBLE PRECISION        :: dr(SATS)             ! rangeの修正量
      DOUBLE PRECISION        :: dx(SATS)             ! 解の更新量
      DOUBLE PRECISION        :: sol(3)               ! 求める方程式の解 receiver positionのx, y, z座標
      INTEGER                 :: i, n, loop, u        ! ループ用カウンタ
      DOUBLE PRECISION        :: x, y, z              ! 解の確認用出力
!   +-----------------------------------------------------------------------------------------------------------------


  SAT1_POS(:) = (/ -13897607.6294d0, -10930188.6233d0, 19676689.6804d0 /)
  SAT2_POS(:) = (/ -17800899.1998d0, 15689920.8120d0, 11943543.3888d0 /)
  SAT3_POS(:) = (/ -1510958.2282d0, 26280096.7818d0, -3117646.1949d0 /)
  SAT4_POS(:) = (/ -12210758.3517d0, 20413597.0201d0, -11649499.5474d0 /)
  SAT5_POS(:) = (/ -170032.6981d0, 17261822.6784d0, 20555984.4061d0 /)


  SATS_POSITION(:, :) = 0.0d0 ! 衛星の位置を初期化
  do i = 1, 3 ! 各衛星の位置を配列にセット
    SATS_POSITION(i, 1) = SAT1_POS(i)
    SATS_POSITION(i, 2) = SAT2_POS(i)
    SATS_POSITION(i, 3) = SAT3_POS(i)
    SATS_POSITION(i, 4) = SAT4_POS(i)
    SATS_POSITION(i, 5) = SAT5_POS(i)
  end do

  ! 観測データを配列にセット
  sats_range(:) = (/ 23634878.5219d0, 20292688.3557d0, 24032055.0372d0, 24383229.3740d0, 22170992.8178d0/)

  sol(:) = 0.d0  ! 解を初期化
  G(:, :) = 0.d0 ! 観測行列を初期化

  ! 解を求めるループ
  do loop=1, MAX_LOOP
    n = SATS  ! 衛星の数をセットする
    do i=1, n ! 衛星の数だけループする
      xyz(1:3) = SATS_POSITION(1:3, i) ! 衛星の位置を配列xyzにセット
      r = sqrt( sum( (xyz(1:3) - sol(1:3) ) ** 2.d0 ) ) ! 疑似距離rの計算
      G(i,1:3) = ( sol(1:3) - xyz(1:3) ) / r ! 観測行列Gを作成

      ! 観測行列のデバッグライト
      ! write(*, *) '========================================='
      ! do u = 1, n
      !   write(*, *) G(u, 1:3)
      ! end do

      dr(i) = SATS_RANGE(i) - r !擬似距離の修正量drを計算
    end do

    call least_squares(G, dr, dx, n, 3); ! 最小二乗法により方程式を解く

    do i=1, 3
      sol(i) = sol(i) + dx(i) ! 解を更新
    end do

    ! 途中経過を出力
    write(6, '("LOOP ",I0, 5X,"x = ",f12.3,5X,"y = ",f12.3,5X,"z = ", f12.3)') &
      loop, sol(1), sol(2), sol(3)

    ! ! 正しい解
    x = -3947762.486d0
    y = 3364401.302d0
    z = 3699431.992d0
  end do

  write(6, *) "*****************************************************"
  write(6, '("x = ",f12.3,5X,"y = ",f12.3,5X,"z = ",f12.3)') x, y, z

end program main
