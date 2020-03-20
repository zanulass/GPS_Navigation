program main
  use compute_solution
  implicit none
  INTEGER, PARAMETER :: SATS = 5
  INTEGER, PARAMETER :: MAX_LOOP = 8
  DOUBLE PRECISION SAT1_POS(3), SAT2_POS(3), SAT3_POS(3), SAT4_POS(3), SAT5_POS(3)
  DOUBLE PRECISION SATS_POSITION(3, SATS), SATS_RANGE(SATS)
  DOUBLE PRECISION g_r, g_xyz(3)
  DOUBLE PRECISION g_G(SATS, 3), g_dr(SATS), g_dx(SATS)
  DOUBLE PRECISION g_sol(3) ! 求める方程式の解 receiver positionのx, y, z座標
  INTEGER i, n, loop, u
  DOUBLE PRECISION ::  x, y, z ! 解の確認用出力

  SAT1_POS(:) = (/ -13897607.6294d0, -10930188.6233d0, 19676689.6804d0 /)
  SAT2_POS(:) = (/ -17800899.1998d0, 15689920.8120d0, 11943543.3888d0 /)
  SAT3_POS(:) = (/ -1510958.2282d0, 26280096.7818d0, -3117646.1949d0 /)
  SAT4_POS(:) = (/ -12210758.3517d0, 20413597.0201d0, -11649499.5474d0 /)
  SAT5_POS(:) = (/ -170032.6981d0, 17261822.6784d0, 20555984.4061d0 /)

  SATS_POSITION(:, :) = 0.0d0
  do i = 1, 3
    SATS_POSITION(i, 1) = SAT1_POS(i)
    SATS_POSITION(i, 2) = SAT2_POS(i)
    SATS_POSITION(i, 3) = SAT3_POS(i)
    SATS_POSITION(i, 4) = SAT4_POS(i)
    SATS_POSITION(i, 5) = SAT5_POS(i)
  end do

  sats_range(:) = (/ 23634878.5219d0, 20292688.3557d0, 24032055.0372d0, 24383229.3740d0, 22170992.8178d0/)

  ! 解を初期化
  g_sol(:) = 0.0d0
  ! 観測行列を初期化
  g_G(:, :) = 0.0d0

  ! 解を求めるループ
  do loop=1, MAX_LOOP
    n = SATS
    do i=1, n
      g_xyz(1:3) = SATS_POSITION(1:3, i)

      ! /* デザイン行列をつくる */
      g_r = sqrt( sum( ( g_xyz(1:3) - g_sol(1:3) )**2 ) )
      g_G(i, 1:3) = ( g_sol(1:3) - g_xyz(1:3) ) / g_r

      ! write(*, *) '============================================='
      ! do u = 1, n
      !   write(*, '(100f12.4)') g_G(u, 1:3)
      ! end do

			! /* 擬似距離の修正量 */
      g_dr(i) = SATS_RANGE(i) - g_r
    end do

		! /* 方程式を解く */
    call least_squares(g_G,g_dr,g_dx,n,3);

		! /* 初期値に加える */
    do i=1, 3
      g_sol(i) = g_sol(i) + g_dx(i)
    end do

		! /* 途中経過を出力する */
    write(6, '("LOOP ",I0, 5X,"x = ",f12.3,5X,"y = ",f12.3,5X,"z = ", f12.3)') &
      loop, g_sol(1), g_sol(2), g_sol(3)

    ! 正しい解
    x = -3947762.486d0
    y = 3364401.302d0
    z = 3699431.992d0
  end do

  write(6, *) "*****************************************************"
  write(6, '("x = ",f12.3,5X,"y = ",f12.3,5X,"z = ",f12.3)') x, y, z

end program main
