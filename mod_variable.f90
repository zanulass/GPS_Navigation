! モジュールの定義
module mod_variable
  implicit none

  ! 定数の定義
  real(8), parameter :: mu = 3.986005d14 ! gravitational constant[m^3/s^2]
  real(8), parameter :: OMEGAe_DOT = 7.2921151467d-5 ! earth's rotation rate [rad/s]

  ! 変数
  ! GPS NAVIGATION MESSAGE FILEのパラメータ
  ! 一行目 PRN / EPOCH / SV CLK
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

end module mod_variable
