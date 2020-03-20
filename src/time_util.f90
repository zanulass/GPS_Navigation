module time_util
  use mod_variable
  implicit none
contains
  subroutine julianday(Y, M, D, hh, mm, ss, MJD_now)
    ! Jurian Dayを計算
    implicit none
    real(8) Y, M, D, hh, mm, ss, MJD, MJD_now
    Y = 2000 + Y

    if (M == 1 .or. M == 2) then
      M = M + 12
      Y = Y - 1
    end if

    MJD = floor(365.25*Y) + floor(Y/400) - floor(Y/100) &
       + floor(30.59*(M-2)) + D - 678912 ! + 1721088.5
    ! write(*, *) " MJD = ", MJD

    MJD_now = MJD + (hh/24 + mm/1440 + ss/86400)
  end subroutine julianday
  subroutine utc_to_GPStime(year, month, day, hour, minute, ss)
    ! 日時からGPS system timeへの変換
    ! GPS_week 週番号, sec_week 週の始めからの経過秒
    implicit none
    integer year, month, day, hour, minute
    real(8) Y, M, D, hh, mm, ss, MJD_now, delta_sec
    MJD_now = 0.0 ! 単位はday
    Y = dble(year)
    M = dble(month)
    D = dble(day)
    hh = dble(hour)
    mm = dble(minute)
    call julianday(Y, M, D, hh, mm, ss, MJD_now)
    delta_sec = MJD_now * 86400 - MJD_GPS_ZERO * 86400 ! 単位は秒
    GPS_week = delta_sec / WEEK_SEC
    GPS_sec = mod(delta_sec, WEEK_SEC)
  end subroutine utc_to_GPStime
end module time_util
