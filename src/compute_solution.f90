module compute_solution
  use mod_variable
  implicit none
contains
  subroutine inverse_matrix(a, m)
  ! ------------------------------------------------------------
  !   a(m, m)   入力される行列,サブルーチン内で逆行列に書き換える
  !   m         行列の次元
  !   b(m, m)   行基本操作を行う操作用の行列
  ! ------------------------------------------------------------
    ! 逆行列を計算する
    implicit none
    INTEGER i, j, k
    INTEGER, INTENT(IN) :: m ! 行列の次元
    DOUBLE PRECISION, INTENT(INOUT) :: a(m, m)
    DOUBLE PRECISION b(MAX_UNKNOWNS, MAX_UNKNOWNS*2)

    ! 操作用の行列をつくる
    do i=1, m
      do j=1, m
        b(i, j) =a(i, j)
        if (i == j) then
          b(i, j+m) =1.0
        else
          b(i, j+m)=0.0
        end if
      end do
    end do

    ! ガウスの消去法
    do i=1, m
  		! 第i行をb[i][i]で正規化する
      if (abs(b(i, i))<=1d-10) stop "Cannot inverse matrix."
      j = m + m
      do while(j >= i)
        b(i, j) = b(i, j) / b(i, i)
        j = j - 1
      end do

  		! 他の行の第i列を消去する
      do k=1, m
        if (k /= i) then
          j = m + m
          do while (j>=i)
            b(k, j) = b(k, j) - (b(k, i) * b(i, j))
            j = j - 1
          end do
        end if
      end do
    end do

  	! 元の行列を逆行列で上書きする
    do i=1, m
      do j=1, m
        a(i, j) = b(i, j + m)
      end do
    end do

    ! write(*, *) '------------------------------------------'
    ! do i = 1, m
    !   write(*, '(100f12.4)') b(i, 1:m)
    ! end do

  end subroutine inverse_matrix

  subroutine least_squares(G, dr, dx, n, m)
  ! ------------------------------------------------------------
  !   G(n, m) 観測行列(n×m)
  !   dr(m)    方程式の右辺(range)
  !   dx(m)    方程式の解で上書きされる(receiver position)
  !   n       方程式の数
  !   m       未知数の数
  ! ------------------------------------------------------------
    integer, intent(in) ::  n, m
    real(8) G(n, m), dr(n), dx(m)
    ! 作業用
    integer i, j
    real(8) GT(m, n), GTG(m, m), GTG_inv(m, m), GTG_invGT(m, n)
    ! 作業用行列を初期化
    GT(:, :) = 0.0d0
    GTG(:, :) = 0.0d0
    GTG_inv(:, :) = 0.0d0
    GTG_invGT(:, :) = 0.0d0

    ! Gの転置行列GT
    GT = transpose(G)

    ! GTとGの行列積
    GTG = matmul(GT, G)

    ! (G^TG)^-1 (逆行列)
    GTG_inv = GTG
    call inverse_matrix(GTG_inv, m)

    ! GTG_invとGの行列積
    GTG_invGT = matmul(GTG_inv, GT)

    ! dx = GTG_invGTdr
    do i = 1, m
      dx(i) = 0.0d0
      do j = 1, n
        dx(i) = dx(i) + GTG_invGT(i, j) * dr(j)
      end do
    end do

  end subroutine least_squares

end module compute_solution
