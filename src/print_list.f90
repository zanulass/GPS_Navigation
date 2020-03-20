program print_list
  use mod_variable
  use file_io
  implicit none

  call read_nav_msg()

  call print_ephemeris_info(ephem_list)

contains
  subroutine print_ephemeris_info(list)
    TYPE(ephemeris_info) list(MAX_EPHMS)
    INTEGER :: i

    do i = 1, size(list)
      ephem_buf = ephem_list(i)
      write(*, *) ephem_buf%PRN, ephem_buf%AF0, ephem_buf%AF1, ephem_buf%AF2
      write(*, *) ephem_buf%TOC
      ! ----------- 2行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%IODE, ephem_buf%Crs, ephem_buf%delta_n, ephem_buf%M0
      ! ----------- 3行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%Cuc, ephem_buf%e,  ephem_buf%Cus, ephem_buf%sqrtA
      ! ----------- 4行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%TOE, ephem_buf%Cic, ephem_buf%LOMEGA0, ephem_buf%Cis
      ! ----------- 5行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%i0, ephem_buf%Crc, ephem_buf%somega, ephem_buf%OMEGA_DOT
      ! ----------- 6行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%IDOT, ephem_buf%CAonL2, ephem_buf%WEEK, ephem_buf%L2P
      ! ----------- 7行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%acc, ephem_buf%health, ephem_buf%TGD, ephem_buf%IODC
      ! ----------- 8行目 読み込み ----------------------------------
      write(6, '(3X, 4D19.12)') ephem_buf%TOT, ephem_buf%Fit
    end do
  end subroutine print_ephemeris_info
end program print_list
