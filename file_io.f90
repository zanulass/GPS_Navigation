module file_io
  use mod_variable
  implicit none
contains
  subroutine read_GPS_Nav()
    ! GPS NAVIGATION MESSAGE FILE からパラメータを読み込む
    ! エフェメリスデータの取り出し番号
    integer, parameter :: &
      EPHM_TOC = 1,	EPHM_AF0 = 2,	EPHM_AF1 = 3,	EPHM_AF2= 4, ! line 1
      EPHM_IODE = 5,	EPHM_Crs = 6,	EPHM_d_n = 7,	EPHM_M0 = 8, ! line 2
      EPHM_Cuc = 9,	EPHM_e = 10, EPHM_Cus = 11,	EPHM_sqrtA = 12, ! line 3
      EPHM_TOE = 13, EPHM_Cic = 14,	EPHM_OMEGA0 = 15, EPHM_Cis = 16, ! line 4
      EPHM_i0 = 17,	EPHM_Crc = 18,	EPHM_omega = 19,	EPHM_dOmega= 20, !line 5
      EPHM_di = 21,	EPHM_CAonL2 = 22, EPHM_WEEK = 23,	EPHM_L2P = 24,	! line 6
      EPHM_acc = 25, EPHM_health = 26, EPHM_TGD = 27,	EPHM_IODC = 28,	! line 7
      EPHM_TOT = 29, EPHM_FIT = 30 ! line 8

    real(8) data(8 * 4) ! １衛星のエフェメリスデータ
    integer i, j, n, prn, line, ios, data_num

    ! ファイルオープン
    open(10, file="data/1222040h.20n")

    ! ヘッダ部 ---------------------------------------
    ! 現時点では読み飛ばす
    !do i = 1, 8 ! RINEX 2.11
    do i = 1, 4 ! RINEX 2.10
      read(10, '()')
    end do
    ! -------------------------------------------------

    ! データ部 ----------------------------------------
    write(6, *) "Reading RINEX Nav. . ."
    do data_num = 1, 3
      ! エフェメリスデータ格納配列を初期化
      ephm_data(8 * 4) = 0.0
      ! １行ずつデータを読み込む
      read(10, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12)') &
        PRN, year, month, day, hour, minute, second, &
        ephm_data(EPHM_AF0), ephm_data(EPHM_AF1), ephm_data(EPHM_AF2)
      GPS_sec = 0.0
      call utc_to_GPStime(year, month, day, hour, minute, second)
      ephm_data(EPHM_TOC) = GPS_sec
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_IODE), ephm_data(EPHM_Crs), ephm_data(EPHM_d_n), &
        ephm_data(EPHM_M0)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_Cuc), ephm_data(EPHM_e), ephm_data(EPHM_Cus), &
        ephm_data(EPHM_sqrtA)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_TOE), ephm_data(EPHM_Cic), ephm_data(EPHM_OMEGA0), &
        ephm_data(EPHM_Cis)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_i0), ephm_data(EPHM_Crc), ephm_data(EPHM_omega), &
        ephm_data(EPHM_OMEGADOT)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_IDOT), ephm_data(EPHM_CAonL2), &
        ephm_data(EPHM_WEEK), ephm_data(EPHM_L2P)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_acc), ephm_data(EPHM_health), ephm_data(EPHM_TGD), &
        ephm_data(EPHM_IODC)
      read(10, '(3X, 4D19.12)') ephm_data(EPHM_TOT), ephm_data(EPHM_FIT), &
        ephm_data(EPHM_spare1), ephm_data(EPHM_spare2)
      ! 1衛星分のエフェメリスデータを格納
      do i = 1, 32
        Nav_data(data_num, i) = ephm_data(i)
      end do
    end do
    close(10)

    ! write(6, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2, F5.1,3D19.12)') &
    !   PRN, year, month, day, hour, minute, second, &
    !   clock_bias, clock_drift, clock_drift_rate
    ! write(6, '(3X, 4D19.12)') IODE, Crs, Delta_n, M0
    ! write(6, '(3X, 4D19.12)') Cuc, e, Cus, sqrt_A
    ! write(6, '(3X, 4D19.12)') TOE, Cic, LOMEGA0, Cis
    ! write(6, '(3X, 4D19.12)') i0, Crc, somega, OMEGA_DOT
    ! write(6, '(3X, 4D19.12)') IDOT, Codes_on_L2_channel, GPS_Week, &
    !   L2_P_data_flag
    ! write(6, '(3X, 4D19.12)') accuracy, health, TGD, IODC_Issue_of_Data
    ! write(6, '(3X, 4D19.12)') Transmission_time_of_message, Fit_interval, &
    !   spare1, spare2

  end subroutine read_GPS_Nav
end module file_io
