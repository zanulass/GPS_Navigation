module matrix
  implicit none
contains
  subroutine inverse_matrix(a, m)
    ! 逆行列を計算する
    implicit none
    integer i, j, k, u
    integer, intent(in) :: m ! 行列の次元
    real(8), intent(inout) :: a(m, m)
    real(8) b(m, m)

    ! write(*, *) '------------------------------------------'
    ! do u = 1, m
    !   write(*, '(100f12.4)') a(u, 1:m)
    ! end do

    ! /* 操作用の行列をつくる */
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

    ! /* ガウスの消去法 */
  	do i=1, m
  		! /* 第i行をb[i][i]で正規化する */
  		if (abs(b(i, i))<=1E-10) stop "Cannot inverse matrix."
      j = m + m
      do while(j >= i)
  			b(i, j) = b(i, j) / b(i, i)
        j = j - 1
  		end do


  		! /* 他の行の第i列を消去する */
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

  	! /* 元の行列を逆行列で上書きする */
  	do i=1, m
  		do j=1, m
  			a(i, j) = b(i, j + m)
  		end do
  	end do
    write(*, *) '------------------------------------------'
    do i = 1, m
      write(*, '(100f12.4)') b(i, 1:m)
    end do
  end subroutine inverse_matrix
end module matrix

! program main
!   use matrix
!   implicit none
!   ! real(8), allocatable :: a(:, :)
!   integer m, i
!   real(8) a(4, 4)
!   a = reshape((/ 3., 2., 1., 7.,  5., 3., 9., 7.,  4., 1., 2., 3.,  4., 5., 6., 7./), (/4, 4/))
!   ! 行列の出力
!   write(*, *) '------------------------------------------'
!   do i = 1, 4
!     write(*, '(100f12.4)') a(i, 1:4)
!   end do
!   m = 4
!   call inverse_matrix(a, m)
!   write(*, *) '------------------------------------------'
!   do i = 1, 4
!     write(*, '(100f12.4)') a(i, 1:4)
!   end do
! end program main
