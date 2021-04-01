module satellite_position
  use WGS84_constants
  use time_util
  use navigation_message
  implicit none
contains
  subroutine calc_sat_position(prn, wt, range, sat_position)
    ! 引数詳細
    INTEGER, INTENT(IN)  :: prn
    TYPE(wtime), INTENT(IN)          :: wt
    DOUBLE PRECISION, INTENT(IN) :: range
    DOUBLE PRECISION, INTENT(OUT) :: sat_position(3)

    ! 使用局所領域
    DOUBLE PRECISION :: tk0 = 0.d0
    DOUBLE PRECISION :: tk = 0.d0
    DOUBLE PRECISION :: A ! Semi-major axis
    DOUBLE PRECISION :: e ! eccentricity
    DOUBLE PRECISION :: i0  ! inclination
    DOUBLE PRECISION :: Ek  ! Eccentric anomaly
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

    !  軌道エポックからの時間差(Time from ephemeris reference epoch)
    tk0 = ( wt%week - current_ephem(prn)%WEEK ) * WEEK_SEC &
          + wt%sec - current_ephem(prn)%TOE
    tk = tk0 - range / C  ! 伝搬時間分を引く

    ! 離心近点角(eccentric anomaly)[rad]を求める
    call calc_eccentric_anomaly(prn, tk, Ek)

    ! 離心率(eccentricity)を取得
    e = current_ephem(prn)%e ! eccentricity

    ! 真近点離角(True anomary)を計算
    vk = atan2( sqrt( 1.d0 - e**2.d0 ) * sin(Ek), cos(Ek) - e )
    ! 緯度引数(Argument of Latitude)[rad]を計算
    PHI_k = vk + current_ephem(prn)%somega !

    ! Second Harmonic Perturbations(補正係数の計算)
    delta_uk = current_ephem(prn)%Cus * sin(2.d0 * PHI_k) &
               + current_ephem(prn)%Cuc * cos(2.d0 * PHI_k) ! Argument of Latitude Correction
    delta_rk = current_ephem(prn)%Crs * sin(2.d0 * PHI_k) &
               + current_ephem(prn)%Crc * cos(2.d0 * PHI_k) ! Radius Correction
    delta_ik = current_ephem(prn)%Cis * sin(2.d0 * PHI_k) &
               + current_ephem(prn)%Cic * cos(2.d0 * PHI_k) ! Inclination Correction

    ! 補正後の緯度引数(Corrected Argument of Latitude)
    uk = PHI_k + delta_uk

    ! 軌道長半径(Semi-major axis)の取得
    A = current_ephem(prn)%sqrtA ** 2
    ! 補正後の軌道長半径(Corrected Radius)
    rk = A * (1.d0 - e * cos(Ek)) + delta_rk

    ! 軌道傾斜角(Inclination)の取得
    i0 = current_ephem(prn)%i0
    ! 補正後の軌道傾斜角(Corrected Inclination)
    ik = i0 + delta_ik + current_ephem(prn)%IDOT * tk

    ! 軌道面内での位置(Positions in orbital plane)
    xk_prime = rk * cos(uk)
    yk_prime = rk * sin(uk)

    ! Corrected longitude of ascending node
    OMEGA_k = current_ephem(prn)%LOMEGA0 + ( current_ephem(prn)%OMEGA_DOT - OMEGAe_DOT ) * tk0 &
              - OMEGAe_DOT * current_ephem(prn)%TOE

    ! Earth-fixed coordinates
    sat_position(1) = xk_prime * cos(OMEGA_k) - yk_prime * cos(ik) * sin(OMEGA_k)
    sat_position(2) = xk_prime * sin(OMEGA_k) + yk_prime * cos(ik) * cos(OMEGA_k)
    sat_position(3) = yk_prime * sin(ik)

  end subroutine calc_sat_position

  subroutine calc_eccentric_anomaly(prn, tk, Ek)
    ! 引数詳細
    INTEGER, INTENT(IN)              :: prn
    DOUBLE PRECISION, INTENT(IN)     :: tk  ! 時間差
    DOUBLE PRECISION, INTENT(INOUT)  :: Ek  ! 離心近点角

    ! 使用局所領域
    DOUBLE PRECISION :: A ! Semi-major axis
    DOUBLE PRECISION :: e ! eccentricity
    DOUBLE PRECISION :: n ! mean motion
    DOUBLE PRECISION :: Mk ! mean anomaly
    DOUBLE PRECISION :: E0, E1  ! Eccentric anomalyのテンポラリ
    DOUBLE PRECISION :: fE_equation ! ケプラー方程式
    DOUBLE PRECISION :: fE_prime ! fEの微分
    INTEGER          :: i  ! ループカウンタ

    ! 離心近点角Ek[rad]を求める
    A = current_ephem(prn)%sqrtA ** 2.d0 ! Semi-major axis
    e = current_ephem(prn)%e ! eccentricity
    n = sqrt(MU / A**3.d0) ! Computed mean motion [rad/sec]
    n = n + current_ephem(prn)%delta_n ! Corrected mean motion
    Mk = current_ephem(prn)%M0 + n * tk ! Mean Anomaly

    ! Solve Kepler equation by Newton-Raphson Method
    E0 = Mk
    ! Ek = Mk
    do i = 1, 10
      ! Ek = Mk + e * sin(Ek)
      fE_equation = E0 - e * sin(E0) - Mk ! 離心近点角Eを変数とする関数f(E)
      fE_prime = 1.d0 - e * cos(E0) ! f(E)の導関数
      E1 = E0 - (fE_equation / fE_prime) ! 反復計算
      if (abs(E1 - E0) < 1d-6) exit ! 許容誤差以内ならexit
      E0 = E1 ! 解の更新
    end do
    Ek = E0 ! 最終的に得られた解(Eccentric Anomaly)

  end subroutine calc_eccentric_anomaly
end module satellite_position
