program solve_kepler_equation
  implicit none
  real(8) E0, Mk, FE, e, dFdE, E1, Ek
  integer i

  Mk = 1.65493151d0
  e = 0.09338096d0

  ! Solve Kepler equation
  E0 = Mk
  do i = 1, 10
    FE = E0 - e * sin(E0) - Mk ! 離心近点角Eを変数とする関数F(E)
    dFdE = 1.0 - e * cos(E0) ! F(E)の導関数

    E1 = E0 - (FE / dFdE) ! 反復計算
    if (abs(E1 - E0) < 1d-6) exit ! 許容誤差以内ならexit
    E0 = E1 ! 解の更新
  end do
  Ek = E0 ! 最終的に得られた解
  write(*, *) Ek

end program solve_kepler_equation
