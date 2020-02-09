module calendar
  implicit none
contains
  subroutine julianday(Y, M, D, h, min, s)
    implicit none
    ! Jurian Dayを計算
    real(8) Y, M, D, h, min, s, MJD

    if (M == 1 .or. M == 2) then
      M = M + 12
      Y = Y - 1
    end if

    MJD = floor(365.25*Y) + floor(Y/400) - floor(Y/100) &
       + floor(30.59*(M-2)) + D - 678912 ! + 1721088.5

    MJD = MJD + (h/24 + min/1440 + s/86400)

  end subroutine julianday

  subroutine GPS_time
    real(8) t ! GPS system time [second]
    ! effective SV PRN code phase time as message transmission time [second]
    real(8) tsv
    real(8) delta_tsv ! SV PRN code pahse time offset [second]

    delta_tsv = af0 + af1


  end subroutine GPS_time
end module calendar
