module satellite_position
  use mod_variable
  implicit none
contains
  subroutine calc_satpos(wt, current_ephem)
    ! 衛星の位置を計算する
    ! 引数詳細
    TYPE(wtime), INTENT(IN)          :: wt
    TYPE(ephemeris_info), INTENT(INOUT) :: current_ephem

    ! 局所領域
    DOUBLE PRECISION :: tk0, tk ! 時間差
    DOUBLE PRECISION :: A ! Semi-major axis
    DOUBLE PRECISION :: e ! eccentricity
    DOUBLE PRECISION :: n ! mean motion
    DOUBLE PRECISION :: Mk ! mean anomaly
    DOUBLE PRECISION :: E0, E1, Ek ! Eccentric anomaly
    DOUBLE PRECISION :: fE ! ケプラー方程式
    DOUBLE PRECISION :: fE_prime ! fEの微分
    DOUBLE PRECISION :: vk ! True anomaly
    DOUBLE PRECISION :: PHI_k ! Argument of Latitude
    DOUBLE PRECISION :: delta_uk ! Argument of Latitude Correction
    DOUBLE PRECISION :: delta_rk ! Radius Correction
    DOUBLE PRECISION :: delta_ik ! Inclination Correction
    DOUBLE PRECISION :: uk ! Corrected Argument of Latitude
    DOUBLE PRECISION :: rk ! Corrected Radius
    DOUBLE PRECISION :: ik ! Corrected Inclination
    DOUBLE PRECISION :: xk_prime, yk_prime ! Positions in orbital plane
    DOUBLE PRECISION :: OMEGA_k ! Corrected longitude of ascending node
    INTEGER          :: i ! ループ用のカウンタ

    !  Time from ephemeris reference epoch
    tk0 = ( wt%week - current_ephem%WEEK ) * WEEK_SEC + wt%sec - current_ephem%TOE
    tk = tk0

    ! 離心近点角Ek[rad]を求める
    A = current_ephem%sqrtA ** 2 ! Semi-major axis
    e = current_ephem%e ! eccentricity
    n = sqrt(MU / A**3) ! Computed mean motion [rad/sec]
    n = n + current_ephem%delta_n ! Corrected mean motion
    Mk = current_ephem%M0 + n * tk0 ! Mean Anomaly

    ! Solve Kepler equation by Newton-Raphson Method
    E0 = Mk
    do i = 1, 10
      fE = E0 - e * sin(E0) - Mk ! 離心近点角Eを変数とする関数f(E)
      fE_prime = 1.d0 - e * cos(E0) ! f(E)の導関数
      E1 = E0 - (fE / fE_prime) ! 反復計算
      if (abs(E1 - E0) < 1d-6) exit ! 許容誤差以内ならexit
      E0 = E1 ! 解の更新
    end do
    Ek = E0 ! 最終的に得られた解(Eccentric Anomaly)

    vk = atan2( sqrt( 1.d0 - e ** 2.d0 ) * sin(Ek), cos(Ek) - 2.d0 ) ! True Anomaly
    PHI_k = vk + current_ephem%somega ! Argument of Latitude

    ! Second Harmonic Perturbations
    delta_uk = current_ephem%Cus * sin(2.d0 * PHI_k) &
               + current_ephem%Cuc * cos(2.d0 * PHI_k) ! Argument of Latitude Correction
    delta_rk = current_ephem%Crs * sin(2 * PHI_k) &
               + current_ephem%Crc * cos(2 * PHI_k) ! Radius Correction
    delta_ik = current_ephem%Cis * sin(2 * PHI_k) &
               + current_ephem%Cic * cos(2 * PHI_k) ! Inclination Correction

    ! 補正係数を適用する
    uk = PHI_k + delta_uk ! Corrected Argument of Latitude
    rk = A * (1.d0 - e * cos(Ek)) + delta_rk ! Corrected Radius
    ik = current_ephem%i0 + delta_ik + current_ephem%IDOT * tk ! Corrected Inclination

    ! Positions in orbital plane(軌道面内での位置)
    xk_prime = rk * cos(uk)
    yk_prime = rk * sin(uk)

    ! Corrected longitude of ascending node
    OMEGA_k = current_ephem%LOMEGA0 + ( current_ephem%OMEGA_DOT - OMEGAe_DOT ) * tk &
              - OMEGAe_DOT * current_ephem%TOE

    ! Earth-fixed coordinates
    current_ephem%pos_xyz(1) = xk_prime * cos(OMEGA_k) - yk_prime * cos(ik) * sin(OMEGA_k)
    current_ephem%pos_xyz(2) = xk_prime * sin(OMEGA_k) + yk_prime * cos(ik) * cos(OMEGA_k)
    current_ephem%pos_xyz(3) = yk_prime * sin(ik)

  end subroutine calc_satpos
end module satellite_position
