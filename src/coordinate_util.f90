module coordinate_util
  use WGS84_constants
  implicit none
contains
  subroutine ecef_to_blh(ecef_pos, blh_pos)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)    :: ecef_pos(3) ! ECEF座標値
    DOUBLE PRECISION, INTENT(INOUT) :: blh_pos(3)  ! 測地座標値(緯度，軽度，高度)

    ! 使用局所領域
    DOUBLE PRECISION  :: b  ! 楕円体の短半径
    DOUBLE PRECISION  :: e  ! 楕円体の離心率
    DOUBLE PRECISION  :: h, p, t, n
    DOUBLE PRECISION  :: lat, lon, height

    ! 測地座標値の初期化
    blh_pos(:) = 0

    ! 原点の場合
    if (ecef_pos(1) == 0 .and. ecef_pos(2) == 0 .and. ecef_pos(3) == 0) then
      goto 8000
    end if

    ! 楕円体のパラメータ
    ! Fe : 扁平率
    ! Re : 長半径
    b = Re * (1.d0 - Fe)  ! 短半径
    e = sqrt(Fe * (2.d0 - Fe))  ! 離心率

    ! 座標変換のためのパラメータ
    h = Re * Re - b * b
    p = sqrt( ecef_pos(1) * ecef_pos(1) + ecef_pos(2) * ecef_pos(2) )
    t = atan2( ecef_pos(3) * Re, p * b )

    ! 測地座標系への変換
    lat = atan2( ecef_pos(3) + h / b * sin(t)**3, p - h / Re * cos(t)**3 )
    n = Re / sqrt( 1.d0 - e**2 * sin(lat)**2 )
    lon = atan2( ecef_pos(2), ecef_pos(1) )
    height = ( p / cos(lat) ) - n

    blh_pos(1) = lat
    blh_pos(2) = lon
    blh_pos(3) = height

    8000 continue

  end subroutine ecef_to_blh

  subroutine ecef_to_enu(ecef_pos, base_pos, enu_pos)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)  :: ecef_pos(3)  ! ECEF座標値
    DOUBLE PRECISION, INTENT(IN)  :: base_pos(3)  ! 基準位置
    DOUBLE PRECISION, INTENT(INOUT)  :: enu_pos(3)  ! ENU座標値

    ! 使用局所領域
    DOUBLE PRECISION  :: blh_base_pos(3)  ! 測地座標値
    DOUBLE PRECISION  :: relative_pos(3)  ! 基準位置からの相対位置
    DOUBLE PRECISION  :: s1, c1, s2, c2

    ! 基準位置からの相対位置を計算
    relative_pos(1:3) = ecef_pos(1:3) - base_pos(1:3)

    ! 基準位置の経緯度を計算
    call ecef_to_blh(base_pos, blh_base_pos)
    s1 = sin(blh_base_pos(2))
    c1 = cos(blh_base_pos(2))
    s2 = sin(blh_base_pos(1))
    c2 = sin(blh_base_pos(1))

    ! 相対位置を回転させてENU座標に変換する
    enu_pos(1) = ecef_pos(1) * s1 * ecef_pos(2) * c1
    enu_pos(2) = ecef_pos(1) * c1 * s2 - ecef_pos(2) * s1 * s2 + ecef_pos(3) * c2
    enu_pos(3) = ecef_pos(3) * c1 * c2 + ecef_pos(2) * s1 * c2 + ecef_pos(3) * s2

  end subroutine ecef_to_enu

  subroutine calc_elevation(sat_position, receiver_position, el)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)  :: sat_position(3)
    DOUBLE PRECISION, INTENT(IN)  :: receiver_position(3)
    DOUBLE PRECISION, INTENT(INOUT)  :: el

    ! 使用局所領域
    DOUBLE PRECISION  :: enu_sat_position(3)

    call ecef_to_enu(sat_position, receiver_position, enu_sat_position)

    el = atan2(enu_sat_position(3), sqrt(enu_sat_position(1)**2 + enu_sat_position(2)**2) )

  end subroutine calc_elevation

  subroutine calc_azimuth(sat_position, receiver_position, az)
    implicit none
    ! 引数詳細
    DOUBLE PRECISION, INTENT(IN)  :: sat_position(3)
    DOUBLE PRECISION, INTENT(IN)  :: receiver_position(3)
    DOUBLE PRECISION, INTENT(INOUT)  :: az

    ! 使用局所領域
    DOUBLE PRECISION  :: enu_sat_position(3)

    call ecef_to_enu(sat_position, receiver_position, enu_sat_position)
    az = atan2(enu_sat_position(1), enu_sat_position(2))

  end subroutine calc_azimuth



end module coordinate_util