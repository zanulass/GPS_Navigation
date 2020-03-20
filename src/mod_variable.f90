　! モジュールの定義
module mod_variable
  implicit none

  ! WGS84定数
  DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535898d0
  DOUBLE PRECISION, PARAMETER :: C = 2.99792458d8
  DOUBLE PRECISION, PARAMETER :: MU = 3.986005d14 ! gravitational constant[m^3/s^2]
  DOUBLE PRECISION, PARAMETER :: OMEGAe_DOT = 7.2921151467d-5 ! earth's rotation rate [rad/s]


  ! GPS NAVIGATION MESSAGE FILEのパラメータ
  ! １行目
  ! integer PRN, year, month, day, hour, minute
  ! DOUBLE PRECISION second, clock_bias, clock_drift, clock_drift_rate
  ! ! 2行目 BROADCAST ORBIT - 1
  ! DOUBLE PRECISION IODE, Crs, Delta_n, M0
  ! ! 3行目 BROADCAST ORBIT - 2
  ! DOUBLE PRECISION Cuc, e, Cus, sqrt_A
  ! ! 4行目 BROADCAST ORBIT - 3
  ! DOUBLE PRECISION TOE, Cic, LOMEGA0, Cis
  ! ! 5行目 BROADCAST ORBIT - 4
  ! DOUBLE PRECISION i0, Crc, somega, OMEGA_DOT
  ! ! 6行目 BROADCAST ORBIT - 5
  ! DOUBLE PRECISION IDOT, Codes_on_L2_channel, GPS_Week, L2_P_data_flag
  ! ! 7行目 BROADCAST ORBIT - 6
  ! DOUBLE PRECISION accuracy, health, TGD, IODC_Issue_of_Data
  ! ! 8行目 BROADCAST ORBIT - 7
  ! DOUBLE PRECISION Transmission_time_of_message, Fit_interval, spare1, spare2


  ! 時刻の定数
  DOUBLE PRECISION, PARAMETER :: WEEK_SEC = 7d0 * 24d0 * 60d0 * 60d0
  INTEGER, PARAMETER :: GPS_ZERO(6) =(/ 1980, 1, 6, 0, 0, 0 /)
  DOUBLE PRECISION, PARAMETER ::  MJD_GPS_ZERO = 44244.00000
  ! GPS timeのパラメータ
  INTEGER GPS_week ! 週番号
  DOUBLE PRECISION GPS_sec  ! 週の始めからの経過秒

  ! 取り扱える行列の大きさ
  INTEGER, PARAMETER :: MAX_SATS = 16 ! 観測衛星数の上限
  INTEGER, PARAMETER :: MAX_UNKNOWNS = 4 ! 未知数の上限
  INTEGER, PARAMETER :: MAX_PRN = 32 ! 衛星番号の上限



end module mod_variable
