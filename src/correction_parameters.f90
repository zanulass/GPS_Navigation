module correction_parameters
    use WGS84_constants
    use time_util
    use coordinate_util
    use navigation_message
    use satellite_position
    implicit none
contains
  subroutine correct_sat_clock(prn, wt, range, sat_clock)
    implicit none
    ! 衛星のクロック誤差を修正
    ! 引数詳細
    INTEGER, INTENT(IN)           :: prn
    TYPE(wtime), INTENT(IN)       :: wt
    DOUBLE PRECISION, INTENT(IN)  :: range
    DOUBLE PRECISION, INTENT(OUT) :: sat_clock

   ! 使用局所領域
    DOUBLE PRECISION  :: tk0 = 0.d0
    DOUBLE PRECISION  :: tk = 0.d0
    DOUBLE PRECISION  :: dt = 0.d0
    DOUBLE PRECISION  :: tr = 0.d0
    DOUBLE PRECISION  :: Ek  ! 離心近点角
    DOUBLE PRECISION  :: e  ! 離心率
    DOUBLE PRECISION  :: sqrtA  ! 軌道傾斜角の平方根

    ! 軌道エポックからの時間差
    tk0 = ( wt%week - current_ephem(prn)%WEEK ) * WEEK_SEC &
          + wt%sec - current_ephem(prn)%TOE
    tk = tk0 - range / C  ! 伝搬時間分を引く

    ! 離心近点角を求める
    call calc_eccentric_anomaly(prn, tk, Ek)
    ! エフェメリスから離心率を取得
    e = current_ephem(prn)%e ! eccentricity
    ! エフェメリスから軌道長半径の平方根を取得
    sqrtA = current_ephem(prn)%sqrtA

    ! 相対論補正項を計算
    tr = ( -2.d0 * sqrt(MU) / C**2 ) * e * sqrtA * sin(Ek)

    ! クロックエポックからの時間差を求める
    tk0 = ( wt%week - current_ephem(prn)%WEEK ) * WEEK_SEC &
          + wt%sec - current_ephem(prn)%TOC
    tk = tk0 - range / C  ! 伝搬時間分を引く

    ! 衛星時計の補正量を計算
    dt = current_ephem(prn)%AF0 &
       + current_ephem(prn)%AF1 * tk &
       + current_ephem(prn)%AF2 * tk * tk

    sat_clock = dt + tr - current_ephem(prn)%TGD

  end subroutine correct_sat_clock

  subroutine calc_iono_correction(sat_position, receiver_position, wt, iono_correction)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)     :: sat_position(3)
    DOUBLE PRECISION, INTENT(IN)     :: receiver_position(3)
    TYPE(wtime), INTENT(IN)          :: wt
    DOUBLE PRECISION, INTENT(INOUT)  :: iono_correction

    ! 使用局所領域
    DOUBLE PRECISION, PARAMETER :: PER_MIN = 20.d0 * 3600.d0
    DOUBLE PRECISION, PARAMETER :: NIGHT_DELAY = 5.d-09
    DOUBLE PRECISION, PARAMETER :: MAX_DELAY_TIME = 14.d0 * 3600.d0
    INTEGER           :: i ! ループカウンタ
    DOUBLE PRECISION, PARAMETER :: RAD_TO_SEC = 1.d0 / PI
    DOUBLE PRECISION  :: el  ! 仰角
    DOUBLE PRECISION  :: az  ! 方位角
    DOUBLE PRECISION  :: blh_receiver_position(3)  ! 測地座標における受信機位置
    DOUBLE PRECISION  :: psi, phi_i, lam_i, phi_m, x, f, y
    DOUBLE PRECISION  :: lt  ! 地方時
    DOUBLE PRECISION  :: amp
    DOUBLE PRECISION  :: per
    DOUBLE PRECISION  :: real_i

    ! 仰角と方位角を求める
    call calc_elevation(sat_position, receiver_position, el)
    call calc_azimuth(sat_position, receiver_position, az)

    ! ピアースポイントを求める
    call ecef_to_blh(receiver_position, blh_receiver_position)
    psi = 0.0137d0 / ( el + 0.11d0 ) - 0.022
    phi_i = RAD_TO_SEC * ( blh_receiver_position(1) + psi * cos(az * PI) )
    if (phi_i > 0.416d0) phi_i = 0.416d0
    if (phi_i < -0.416d0) phi_i = -0.416d0
    lam_i = RAD_TO_SEC * ( blh_receiver_position(2) + psi * sin(az * PI) / cos(phi_i * PI) )
    phi_m = phi_i + 0.064d0 * cos( (lam_i - 1.617d0) * PI )

    ! 地方時
    lt = DAY_SEC / 2.d0 * lam_i + wt%sec
    do while (lt > DAY_SEC)
      lt = lt - DAY_SEC
    end do
    do while (lt < 0.d0)
      lt = lt + DAY_SEC
    end do

    ! コサイン関数の周期と振幅
    amp = 0.d0
    per = 0.d0
    do i=1, IONO_PARAMETERS
      real_i = real(i, kind=8)
      amp = amp + ( ion_alpha(i) * phi_m ** real_i )
      per = per + ( ion_beta(i) * phi_m ** real_i )
    end do

    if (amp < 0.d0) amp = 0.d0
    if (per < PER_MIN) per = PER_MIN

    ! 傾斜係数
    x = 0.53d0 - el
    f = 1.d0 + 16.d0 * x**3

    ! 補正値を求める
    x = 2.d0 * PI * ( lt - MAX_DELAY_TIME ) / per
    do while (x > PI)
      x = x - 2.d0 * PI
    end do
    do while (x < -PI)
      x = x + 2.d0 * PI
    end do
    if (abs(x) < 1.57d0) then
      y = amp * ( 1.d0 - x**2 * (0.5d0 - x**2 / 24.d0) ) ! 昼間
    else
      y = 0.d0  ! 夜間
    end if

    iono_correction = -(f * (NIGHT_DELAY + y) * C) ! 距離にして返す

  end subroutine calc_iono_correction

  subroutine calc_tropo_correction(sat_position, receiver_position, tropo_correction)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)  :: sat_position(3)
    DOUBLE PRECISION, INTENT(IN)  :: receiver_position(3)
    DOUBLE PRECISION, INTENT(INOUT) :: tropo_correction

    ! 使用局所領域
    DOUBLE PRECISION, PARAMETER :: TROPO_DELAY_ZENITH = 2.47d0
    DOUBLE PRECISION, PARAMETER :: TROPO_SCALE_HEIGHT = 1.d0 / 2.3d-5
    DOUBLE PRECISION  :: blh_receiver_position(3)  ! 測地座標における受信機位置
    DOUBLE PRECISION  :: d
    DOUBLE PRECISION  :: el

    call ecef_to_blh(receiver_position, blh_receiver_position)

    if ( blh_receiver_position(3) < 0.d0) then
      d = 1.d0
    else if (blh_receiver_position(3) > TROPO_SCALE_HEIGHT) then
      d = 0.d0
    else
      d = 1.d0 - blh_receiver_position(3) / TROPO_SCALE_HEIGHT
    end if

    ! 仰角を求める
    call calc_elevation(sat_position, receiver_position, el)

    tropo_correction = -(TROPO_DELAY_ZENITH * d**5 / ( sin(el) + 0.0121d0 ) )

  end subroutine

end module
