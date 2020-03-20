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

    write(*, *) '##### input Navigation Message #####'
    do i = 1, size(list)
      ephem_buf = ephem_list(i)
      if (ephem_buf%PRN == 0) exit

      write(*, *) ''
      write(*, '(A,1X,I2)') 'PRN:', ephem_buf%PRN
      ! ----------- 1行目 ----------------------------------
      write(*, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'TOC:', ephem_buf%TOC, 'AF0:', ephem_buf%AF0, 'AF1:', ephem_buf%AF1, &
        'AF2:', ephem_buf%AF2
      ! ----------- 2行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'IODE:', ephem_buf%IODE, 'Crs:', ephem_buf%Crs, 'delta_n:', ephem_buf%delta_n, &
        'M0: ', ephem_buf%M0
      ! ----------- 3行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'Cuc:', ephem_buf%Cuc, 'e:', ephem_buf%e,  'Cus:', ephem_buf%Cus, &
        'sqrtA:', ephem_buf%sqrtA
      ! ----------- 4行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'TOE:', ephem_buf%TOE, 'Cic:', ephem_buf%Cic, 'OMEGA0:', ephem_buf%LOMEGA0, &
        'Cis:', ephem_buf%Cis
      ! ----------- 5行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'i0:', ephem_buf%i0, 'Crc:', ephem_buf%Crc, 'omega:', ephem_buf%somega, &
       'OMEGADOT:', ephem_buf%OMEGA_DOT
      ! ----------- 6行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'IDOT:', ephem_buf%IDOT, 'CAonL2:', ephem_buf%CAonL2, 'WEEK:', ephem_buf%WEEK, &
       'L2P:', ephem_buf%L2P
      ! ----------- 7行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'acc:', ephem_buf%acc, 'health:', ephem_buf%health, 'TGD:', ephem_buf%TGD, &
       'IODC:', ephem_buf%IODC
      ! ----------- 8行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'TOT:', ephem_buf%TOT, 'Fit:', ephem_buf%Fit
      write(*, *) ''
    end do
  end subroutine print_ephemeris_info
end program print_list
