module compute_solution
  use matrix
  implicit none
contains
  subroutine least_squares(G, dr, dx, n, m)
  ! ------------------------------------------------------------
  !   G[][];   デザイン行列(n×m)
  !   dr[];    方程式の右辺(m次)
  !   dx[];    方程式の解で上書きされる(m次)
  !   n;       方程式の数
  !   m;       未知数の数
  ! ------------------------------------------------------------
  integer, intent(in) ::  n, m
  real(8) G(n, m), dr(n), dx(m)
  ! 作業用
  integer i, j, k, u
  real(8) GT(m, n), GTG(m, m), GTG_inv(m, m), GTG_invGT(m, n)

    ! 最小二乗法 dx = (G^T*G)^-1 * G^T * dr

    ! Gの転置行列 G^T
    GT(m, n) = 0.0
    do i = 1, m
      do j = 1, n
        GT(i, j) = G(j, i)
      end do
    end do

    ! G^TG (行列積)
    GTG(m, m) = 0.0
    do i = 1, m
      do j = 1, m
        do k = 1, n
          GTG(i, j) = GT(i, k) * G(k, j)
        end do
      end do
    end do

    ! (G^TG)^-1 (逆行列)
    GTG_inv = GTG


    call inverse_matrix(GTG_inv, m)

    ! GTG_invG^T (行列積)
    GTG_invGT(m, n) = 0.0
    do i = 1, m
      do j = 1, n
        do k = 1, m
           GTG_invGT(i, j) = GTG_inv(i, k) * GT(k, j)
        end do
      end do
    end do

    ! dx = GTG_invGTdr
    do i = 1, m
      dx(i) = 0.0
      do j = 1, n
        dx(i) = dx(i) + GTG_invGT(i, j) * dr(j)
      end do
    end do

  end subroutine least_squares

end module compute_solution
