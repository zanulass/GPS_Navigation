program main
  use mod_variable
  use file_io
  use calc_satpos
  implicit none
  real(8) t ! GPS system tim

  call read_GPS_Nav()
  t = 0
  call ephemeris_determination(t)

end program main
