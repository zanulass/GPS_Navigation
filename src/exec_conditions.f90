module exec_conditions
  ! 扱える行列の大きさ
  INTEGER, PARAMETER :: MAX_SATS = 5 ! 観測衛星数の上限
  INTEGER, PARAMETER :: MAX_UNKNOWNS = 4 ! 未知数の上限
  INTEGER, PARAMETER :: MAX_PRN = 32 ! 衛星番号の上限
  INTEGER, PARAMETER :: MAX_EPHMS = 100 ! 記録するエフェメリスの上限
  DOUBLE PRECISION, PARAMETER :: EPHEMERIS_EXPIRE = 2.d0 ! エフェメリスの有効期限 [h]
  INTEGER, PARAMETER :: MAX_LOOP = 8

  ! ファイルパス
  ! tmp領域
  CHARACTER(256)  :: tmp_dir = "../tmp/"
  ! RINEX NAVIGATION MESSAGE FILEのパス
  CHARACTER(256) :: nav_msg_file = "../data/mtka3180.05n"
  ! 実行結果リストファイルのパス
  CHARACTER(256) :: list_file = "../tmp/list"

  ! 電離層遅延補正フラグ
  LOGICAL  :: iono_flag = .true.
  ! 対流圏遅延補正フラグ
  LOGICAL  :: tropo_flag = .false.

  ! 測位計算に使用した衛星のPRN
  INTEGER  :: num_used_PRN
  INTEGER  :: used_PRN_list(MAX_PRN)

  ! 観測データ補正値記録用csvデータ
  DOUBLE PRECISION  :: sat_clock_for_print(MAX_PRN, MAX_LOOP)
  DOUBLE PRECISION  :: iono_correction_for_print(MAX_PRN, MAX_LOOP)
  DOUBLE PRECISION  :: tropo_correction_for_print(MAX_PRN, MAX_LOOP)


end module exec_conditions