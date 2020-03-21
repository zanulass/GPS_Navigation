module print_list
  use mod_variable
  use navigation_message
  implicit none
  CHARACTER(256)          :: nav_msg_file
  INTEGER, PARAMETER      :: PRN_list(5) = (/ 5,14,16,22,25 /)

  nav_msg_file = "../data/mtka3180.05n"
  call read_nav_msg(nav_msg_file)

  call print_ephemeris_info(PRN_list)

contains
  subroutine print_ephemeris_info(PRN_li)
    INTEGER, INTENT(IN) :: PRN_li(5)
    INTEGER :: i
    TYPE(wtime)             :: wt                   ! 時刻
    INTEGER                 :: iode
    TYPE(ephemeris_info) :: currnt_ephem ! 作業用のエフェメリス

    wt%week = 1349     ! 05/11/13〜19の週
    wt%sec = 86400.d0  ! 月曜日の00:00:00
    iode = -1

    write(*, *) ''
    write(*, *) '##### input Navigation Message #####'
    write(*, *) ''
    write(*, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'alpha1:', ion_alpha(1), 'alpha2:', ion_alpha(2), 'alpha3:', ion_alpha(3), 'alpha4:', ion_alpha(4)
    write(*, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'beta1:', ion_beta(1), 'beta2:', ion_beta(2), 'beta3:', ion_beta(3), 'beta4:', ion_beta(4)
    write(*, '(3X,A9,I6)') 'leap_sec:', leap_sec

    do i = 1, size(PRN_li)
      call set_ephemeris(PRN_li(i), wt, iode, currnt_ephem)
      ! if (currnt_ephem%PRN == 0) exit

      write(*, *) ''
      write(*, '(A,1X,I2)') 'PRN:', currnt_ephem%PRN
      ! ----------- 1行目 ----------------------------------
      write(*, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'TOC:', currnt_ephem%TOC, 'AF0:', currnt_ephem%AF0, 'AF1:', currnt_ephem%AF1, &
        'AF2:', currnt_ephem%AF2
      ! ----------- 2行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'IODE:', currnt_ephem%IODE, 'Crs:', currnt_ephem%Crs, 'delta_n:', currnt_ephem%delta_n, &
        'M0: ', currnt_ephem%M0
      ! ----------- 3行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'Cuc:', currnt_ephem%Cuc, 'e:', currnt_ephem%e,  'Cus:', currnt_ephem%Cus, &
        'sqrtA:', currnt_ephem%sqrtA
      ! ----------- 4行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
        'TOE:', currnt_ephem%TOE, 'Cic:', currnt_ephem%Cic, 'OMEGA0:', currnt_ephem%LOMEGA0, &
        'Cis:', currnt_ephem%Cis
      ! ----------- 5行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'i0:', currnt_ephem%i0, 'Crc:', currnt_ephem%Crc, 'omega:', currnt_ephem%somega, &
       'OMEGADOT:', currnt_ephem%OMEGA_DOT
      ! ----------- 6行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'IDOT:', currnt_ephem%IDOT, 'CAonL2:', currnt_ephem%CAonL2, 'WEEK:', currnt_ephem%WEEK, &
       'L2P:', currnt_ephem%L2P
      ! ----------- 7行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'acc:', currnt_ephem%acc, 'health:', currnt_ephem%health, 'TGD:', currnt_ephem%TGD, &
       'IODC:', currnt_ephem%IODC
      ! ----------- 8行目 ----------------------------------
      write(*, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
       'TOT:', currnt_ephem%TOT, 'Fit:', currnt_ephem%Fit
      write(*, *) ''
    end do
  end subroutine print_ephemeris_info
end module print_list
