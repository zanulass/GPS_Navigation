module time_util
  use mod_variable
  implicit none
contains
  subroutine julianday(Y, M, D, hh, mm, ss, MJD_now)
    ! Jurian Dayを計算
    implicit none
    DOUBLE PRECISION Y, M, D, hh, mm, ss, MJD, MJD_now
    Y = 2000.d0 + Y

    if (M == 1.d0 .or. M == 2.d0) then
      M = M + 12.d0
      Y = Y - 1.d0
    end if

    MJD = floor(365.25d0*Y) + floor(Y/400.d0) - floor(Y/100.d0) &
       + floor(30.59d0*(M-2.d0)) + D - 678912.d0 ! + 1721088.5
    ! write(*, *) " MJD = ", MJD

    MJD_now = MJD + (hh/24.d0 + mm/1440.d0 + ss/86400.d0)
  end subroutine julianday


  subroutine utc_to_GPStime(year, month, day, hour, minute, ss)
    ! 日時からGPS system timeへの変換
    ! GPS_week 週番号, sec_week 週の始めからの経過秒
    implicit none
    INTEGER year, month, day, hour, minute
    DOUBLE PRECISION Y, M, D, hh, mm, ss, MJD_now, delta_sec
    MJD_now = 0.0 ! 単位はday
    Y = dble(year)
    M = dble(month)
    D = dble(day)
    hh = dble(hour)
    mm = dble(minute)
    call julianday(Y, M, D, hh, mm, ss, MJD_now)
    delta_sec = MJD_now * 86400.d0 - MJD_GPS_ZERO * 86400.d0 ! 単位は秒
    GPS_week = int(delta_sec / WEEK_SEC)
    GPS_sec = mod(delta_sec, WEEK_SEC)
  end subroutine utc_to_GPStime
end module time_util
