module print_list
  use mod_variable
  implicit none
contains
  subroutine print_nav_file_header()
    write(20, *) ''
    write(20, *) '##### input Navigation Message #####'
    write(20, *) ''
    write(20, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'alpha1:', ion_alpha(1), 'alpha2:', ion_alpha(2), 'alpha3:', ion_alpha(3), 'alpha4:', ion_alpha(4)
    write(20, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'beta1:', ion_beta(1), 'beta2:', ion_beta(2), 'beta3:', ion_beta(3), 'beta4:', ion_beta(4)
    write(20, '(3X,A9,I6)') 'leap_sec:', leap_sec

  end subroutine print_nav_file_header

  subroutine print_ephemeris_info(current_ephem)
    TYPE(ephemeris_info), INTENT(IN) :: current_ephem

    write(20, *) ''
    write(20, '(A,1X,I2)') 'PRN:', current_ephem%PRN
    ! ----------- 1行目 ----------------------------------
    write(20, '(3X,A9,1X,F19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'TOC:', current_ephem%TOC, 'AF0:', current_ephem%AF0, 'AF1:', current_ephem%AF1, &
      'AF2:', current_ephem%AF2
    ! ----------- 2行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'IODE:', current_ephem%IODE, 'Crs:', current_ephem%Crs, 'delta_n:', current_ephem%delta_n, &
      'M0: ', current_ephem%M0
    ! ----------- 3行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'Cuc:', current_ephem%Cuc, 'e:', current_ephem%e,  'Cus:', current_ephem%Cus, &
      'sqrtA:', current_ephem%sqrtA
    ! ----------- 4行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'TOE:', current_ephem%TOE, 'Cic:', current_ephem%Cic, 'OMEGA0:', current_ephem%LOMEGA0, &
      'Cis:', current_ephem%Cis
    ! ----------- 5行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'i0:', current_ephem%i0, 'Crc:', current_ephem%Crc, 'omega:', current_ephem%somega, &
      'OMEGADOT:', current_ephem%OMEGA_DOT
    ! ----------- 6行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'IDOT:', current_ephem%IDOT, 'CAonL2:', current_ephem%CAonL2, 'WEEK:', current_ephem%WEEK, &
      'L2P:', current_ephem%L2P
    ! ----------- 7行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'acc:', current_ephem%acc, 'health:', current_ephem%health, 'TGD:', current_ephem%TGD, &
      'IODC:', current_ephem%IODC
    ! ----------- 8行目 ----------------------------------
    write(20, '(3X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'TOT:', current_ephem%TOT, 'Fit:', current_ephem%Fit
    write(20, *) ''
    write(20, '(3X,A29,2X,A9,1X,D19.12,2X,A9,1X,D19.12,2X,A9,1X,D19.12)') &
      'Satellite Positon:', 'x:', current_ephem%pos_xyz(1), 'y:', current_ephem%pos_xyz(2), 'z:', current_ephem%pos_xyz(3)

  end subroutine print_ephemeris_info
end module print_list
