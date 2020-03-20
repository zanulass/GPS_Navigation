module file_io
  use mod_variable
  use time_util
  implicit none
contains
  subroutine read_nav_msg()
    ! GPS NAVIGATION MESSAGE FILE からパラメータを読み込む

    CHARACTER(256) :: nav_msg_file
    INTEGER        :: ios ! ファイル読み込みステータス
    INTEGER :: year, month, day, hour, minute ! TOC計算用
    DOUBLE PRECISION :: second ! TOC計算用
    INTEGER          :: i ! ループ用カウンタ




    ! ファイルオープン
    nav_msg_file = "../data/mtka3180.05n"
    open(10, file=nav_msg_file)

    write(*, *) 'Reading RINEX Nav...'

    ! ヘッダ部 ---------------------------------------
    read(10, '()') ! 1行目は読み飛ばす
    read(10, '()') ! 2行目は読み飛ばす
    read(10, '(2X,4D12.4)') ion_alpha(1), ion_alpha(2), ion_alpha(3), ion_alpha(4)
    read(10, '(2X,4D12.4)') ion_beta(1), ion_beta(2), ion_beta(3), ion_beta(4)
    read(10, '()') ! DELTA-UTCは読み飛ばす
    read(10, '(I6)') leap_sec
    read(10, '()') ! ヘッダ最終行
    ! -------------------------------------------------

    ! データ部 読み込み----------------------------------------
    ios = 1 ! ファイル読み込みのiostatを初期化

    do i = 1, 1000 ! ファイルの最終行まで読み込む
      ! ----------- 1行目 読み込み ----------------------------------
      read(10, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12)', iostat = ios) &
        ephem_buf%PRN, year, month, day, hour, minute, second, &
        ephem_buf%AF0, ephem_buf%AF1, ephem_buf%AF2

      if (ios < 0) exit ! ファイル最終行に来たら読み込み終了

      GPS_sec = 0.0 ! 週の始めから経過秒を初期化
      call utc_to_GPStime(year, month, day, hour, minute, second) ! TOCを計算する
      ephem_buf%TOC = GPS_sec

      ! ----------- 2行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%IODE, ephem_buf%Crs, ephem_buf%delta_n, ephem_buf%M0
      ! ----------- 3行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%Cuc, ephem_buf%e,  ephem_buf%Cus, ephem_buf%sqrtA
      ! ----------- 4行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%TOE, ephem_buf%Cic, ephem_buf%LOMEGA0, ephem_buf%Cis
      ! ----------- 5行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%i0, ephem_buf%Crc, ephem_buf%somega, ephem_buf%OMEGA_DOT
      ! ----------- 6行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%IDOT, ephem_buf%CAonL2, ephem_buf%WEEK, ephem_buf%L2P
      ! ----------- 7行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_buf%acc, ephem_buf%health, ephem_buf%TGD, ephem_buf%IODC
      ! ----------- 8行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)', iostat=ios) ephem_buf%TOT, ephem_buf%Fit

      ephem_list(i) = ephem_buf ! エフェメリスを配列に格納

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

  end subroutine read_nav_msg
end module file_io
