module calc_satpos
  use mod_variable
  use calendar
  implicit none
contains
  subroutine ephemeris_determination()
    A = sqrt_A**2 ! Semi-major axis
    n0 = sqrt(mu / A**3) ! Computed mean motion [rad/sec]

    ! epoch計算
    t =
    tk =


  end subroutine ephemeris_determination
end module calc_satpos
