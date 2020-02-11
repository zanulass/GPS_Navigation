! モジュールの定義
module mod_variable
  implicit none

  ! WGS84定数
  real(8), parameter :: PI = 3.1415926535898d0
  real(9), parameter :: C = 2.99792458d8
  real(8), parameter :: MU = 3.986005d14 ! gravitational constant[m^3/s^2]
  real(8), parameter :: OMEGAe_DOT = 7.2921151467d-5 ! earth's rotation rate [rad/s]





  ! GPS NAVIGATION MESSAGE FILEのパラメータ
  ! １行目
  integer PRN, year, month, day, hour, minute
  real(8) second, clock_bias, clock_drift, clock_drift_rate
  ! 2行目 BROADCAST ORBIT - 1
  real(8) IODE, Crs, Delta_n, M0
  ! 3行目 BROADCAST ORBIT - 2
  real(8) Cuc, e, Cus, sqrt_A
  ! 4行目 BROADCAST ORBIT - 3
  real(8) TOE, Cic, LOMEGA0, Cis
  ! 5行目 BROADCAST ORBIT - 4
  real(8) i0, Crc, somega, OMEGA_DOT
  ! 6行目 BROADCAST ORBIT - 5
  real(8) IDOT, Codes_on_L2_channel, GPS_Week, L2_P_data_flag
  ! 7行目 BROADCAST ORBIT - 6
  real(8) accuracy, health, TGD, IODC_Issue_of_Data
  ! 8行目 BROADCAST ORBIT - 7
  real(8) Transmission_time_of_message, Fit_interval, spare1, spare2

  ! GPS timeのパラメータ
  ! 時刻の定数
  real(8) WEEK_SEC = 7d0 * 24d0 * 60d0 * 60d0
  integer GPS_ZERO(6) =(/ 1980, 1, 6, 0, 0, 0 /)
  real(8) MJD_GPS_ZERO = 44244.00000
  integer GPS_week ! 週番号
  real(8) GPS_sec  ! 週の始めからの経過秒


end module mod_variable
