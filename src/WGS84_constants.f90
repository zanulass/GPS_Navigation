module WGS84_constants
  implicit none

  ! WGS84定数
  DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535898d0
  DOUBLE PRECISION, PARAMETER :: C = 2.99792458d8
  DOUBLE PRECISION, PARAMETER :: MU = 3.986005d14 ! gravitational constant[m^3/s^2]
  DOUBLE PRECISION, PARAMETER :: OMEGAe_DOT = 7.2921151467d-5 ! earth's rotation rate [rad/s]
  DOUBLE PRECISION, PARAMETER :: Re = 6378137.d0  ! 地球半径[m]
  DOUBLE PRECISION, PARAMETER :: Fe = 1.d0/298.257223563d0 ! 地球の扁平率
end module WGS84_constants
