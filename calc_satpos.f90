module calc_satpos
  use mod_variable
  use calendar
  implicit none
contains
  subroutine ephemeris_determination()
    A = sqrt_A**2 ! Semi-major axis
    n0 = sqrt(mu / A**3) ! Computed mean motion [rad/sec]

    tk = t - TOE ! Time from ephemeris reference epoch
    if (tk > 302400) then
      tk = tk - 604800
    else if (tk < -302400) then
      tk = tk + 604800

    n = n0 + Delta_n ! Corrected mean motion
    Mk = M0 + n * tk ! Mean Anomaly

    ! Solve Kepler equation by Newton-Raphson Method
    E0 = Mk
    do i = 1, 10
      FE = E0 - e * sin(E0) - Mk ! 離心近点角Eを変数とする関数F(E)
      dFdE = 1.0 - e * cos(E0) ! F(E)の導関数
      E1 = E0 - (FE / dFdE) ! 反復計算
      if (abs(E1 - E0) < 1d-6) exit ! 許容誤差以内ならexit
      E0 = E1 ! 解の更新
    end do
    Ek = E0 ! 最終的に得られた解(Eccentric Anomaly)

    sinvk = sqrt(1 - e**2) * sin(Ek) / (1 - e * cos(Ek))
    cosvk = cos(Ek - e) / (1 -e * cos(Ek))
    vk = atan(sinvk / cosvk) ! True Anomaly
    PHI_k = vk + somega ! Argument of Latitude

    ! Second Harmonic Perturbations
    ! Argument of Latitude Correction
    delta_uk = Cus * sin(2 * PHI_k) + Cuc * cos(2 * PHI_k)
    ! Radius Correction
    delta_rk = Crs * sin(2 * PHI_k) + Crc * cos(2 * PHI_k)
    ! Inclination Correction
    delta_ik = Cis * sin(2 * PHI_k) + Cic * cos(2 * PHI_k)



  end subroutine ephemeris_determination
end module calc_satpos
