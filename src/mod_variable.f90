module mod_variable
  implicit none

  ! WGS84定数
  DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535898d0
  DOUBLE PRECISION, PARAMETER :: C = 2.99792458d8
  DOUBLE PRECISION, PARAMETER :: MU = 3.986005d14 ! gravitational constant[m^3/s^2]
  DOUBLE PRECISION, PARAMETER :: OMEGAe_DOT = 7.2921151467d-5 ! earth's rotation rate [rad/s]

  ! 時刻の定数
  DOUBLE PRECISION, PARAMETER :: WEEK_SEC = 7.d0 * 24.d0 * 60.d0 * 60.d0
  INTEGER, PARAMETER :: GPS_ZERO(6) =(/ 1980, 1, 6, 0, 0, 0 /)
  DOUBLE PRECISION, PARAMETER ::  MJD_GPS_ZERO = 44244.00000

  ! 時刻を表す構造体
  type :: wtime
    INTEGER          :: week ! 週番号
    DOUBLE PRECISION :: sec  ! 週初めからの経過時間
  end type wtime

  INTEGER, PARAMETER :: MAX_SATS = 16 ! 観測衛星数の上限
  INTEGER, PARAMETER :: MAX_UNKNOWNS = 4 ! 未知数の上限
  INTEGER, PARAMETER :: MAX_PRN = 32 ! 衛星番号の上限
  INTEGER, PARAMETER :: MAX_EPHMS = 100 ! 記録するエフェメリスの上限
  DOUBLE PRECISION, PARAMETER :: EPHEMERIS_EXPIRE = 2.d0 ! エフェメリスの有効期限 [h]

  ! Navigation Message File ヘッダ部
  DOUBLE PRECISION :: ion_alpha(4)
  DOUBLE PRECISION :: ion_beta(4)
  INTEGER          :: leap_sec
  type :: ephemeris_info
  !----------- 1行目 ------------------------
    INTEGER           :: PRN = 0
    DOUBLE PRECISION  :: TOC = 0.d0, AF0 = 0.d0, AF1 = 0.d0, AF2 = 0.d0
  !----------- 2行目 -------------------------
    DOUBLE PRECISION  :: IODE = 0.d0, Crs = 0.d0, delta_n = 0.d0, M0 = 0.d0
  !----------- 3行目 -------------------------
    DOUBLE PRECISION  :: Cuc = 0.d0, e = 0.d0, Cus = 0.d0, sqrtA = 0.d0
  !----------- 4行目 -------------------------
    DOUBLE PRECISION  :: TOE = 0.d0, Cic = 0.d0, LOMEGA0 = 0.d0, Cis =0.d0
  !----------- 5行目 -------------------------
    DOUBLE PRECISION  :: i0 = 0.d0, Crc = 0.d0, somega = 0.d0, OMEGA_DOT =0.d0
  !----------- 6行目 -------------------------
    DOUBLE PRECISION  :: IDOT = 0.d0, CAonL2 = 0.d0, WEEK = 0.d0, L2P = 0.d0
  !----------- 7行目 -------------------------
    DOUBLE PRECISION  :: acc = 0.d0, health = 0.d0, TGD = 0.d0, IODC = 0.d0
  !----------- 8行目 -------------------------
    DOUBLE PRECISION  :: TOT = 0.d0, Fit = 0.d0
  !----------- ECEF座標における衛星位置
    DOUBLE PRECISION  :: pos_xyz(3)
  end type ephemeris_info

  TYPE(ephemeris_info) :: ephem_buf(MAX_PRN, MAX_EPHMS) ! 全エフェメリス格納配列
  TYPE(ephemeris_info) :: ephem_data ! 1衛星分のエフェメリス
  INTEGER              :: ephem_count(MAX_PRN) = 0




end module mod_variable
