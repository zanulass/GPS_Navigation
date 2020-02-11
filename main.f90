program main
  use matrix
  use compute_solution
  implicit none
  integer, parameter :: SATS = 5
  integer, parameter :: MAX_LOOP = 8
  real(8) sat1_pos(3), sat2_pos(3), sat3_pos(3), sat4_pos(3), sat5_pos(3)
  real(8) sat1_range, sat2_range, sat3_range, sat4_range, sat5_range
  real(8) sats_position(5, 3), sats_range(5)
  integer i, n, loop, u
  real(8) r, x, y, z
  real(8) G(4, 4), dr(4), dx(4)
  real(8) sol(4)

  sat1_pos = (/ -13897607.6294,-10930188.6233,19676689.6804 /)
  sat2_pos = (/ -17800899.1998,15689920.8120,11943543.3888 /)
  sat3_pos = (/ -1510958.2282,26280096.7818,-3117646.1949 /)
  sat4_pos = (/ -12210758.3517,20413597.0201,-11649499.5474 /)
  sat5_pos = (/ -170032.6981,17261822.6784,20555984.4061 /)

  sats_position(5, 3) = 0.0
  do i = 1, 3
    sats_position(1, i) = sat1_pos(i)
    sats_position(2, i) = sat2_pos(i)
    sats_position(3, i) = sat3_pos(i)
    sats_position(4, i) = sat4_pos(i)
    sats_position(5, i) = sat5_pos(i)
  end do

  sats_range = (/ 23634878.5219, 20292688.3557, 24032055.0372, 24383229.3740, 22170992.8178/)

  ! 解を初期化
  sol(3) = 0.0

  ! 解を求めるループ
  do loop = 1, MAX_LOOP
    n = SATS ! 衛星の数
    do i = 1, n
      x = sats_position(i, 1)
      y = sats_position(i, 2)
      z = sats_position(i, 3)
      ! デザイン行列をつくる
      r = sqrt((x - sol(1))**2 + (y - sol(2))**2 + (z - sol(3))**2)
      G(i, 1) = (sol(1) - x) / r
      G(i, 2) = (sol(2) - y) / r
      G(i, 3) = (sol(3) - z) / r

      write(*, *) '------------------------------------------'
      do u = 1, 5
        write(*, '(100f12.4)') G(u, 1:3)
      end do

      ! 疑似距離の修正量 o - c (c = r + bias(今はbiasは考慮しない))
      dr(i) = sats_range(i) - r
    end do

    ! 方程式を解く
    call least_squares(G, dr, dx, n, 3)

    ! 初期値に加える
    do i = 1, 3
      sol(i) = sol(i) + dx(i)
    end do

    ! 途中経過を出力する
    write(*, *) "LOOP = ", loop, "x = ", sol(1), "y = ", sol(2), "z = ", sol(3)
  end do


end program main
