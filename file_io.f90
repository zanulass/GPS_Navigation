module file_io
  use mod_variable
  implicit none
contains
  subroutine read_GPS_Nav()
    ! GPS NAVIGATION MESSAGE FILE からパラメータを読み込む

    integer i, ios ! ループ用カウンタ iostat用変数

    open(10, file="data/brdc0010.19n")
    ! ヘッダ部を読みとばす
    do i=1, 8
      read(10, '()')
    end do
    read(10, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12)') &
      PRN, year, month, day, hour, minute, second, &
      clock_bias, clock_drift, clock_drift_rate
    read(10, '(3X, 4D19.12)') IODE, Crs, Delta_n, M0
    read(10, '(3X, 4D19.12)') Cuc, e, Cus, sqrt_A
    read(10, '(3X, 4D19.12)') TOE, Cic, LOMEGA0, Cis
    read(10, '(3X, 4D19.12)') i0, Crc, somega, OMEGA_DOT
    read(10, '(3X, 4D19.12)') IDOT, Codes_on_L2_channel, GPS_Week, &
      L2_P_data_flag
    read(10, '(3X, 4D19.12)') accuracy, health, TGD, IODC_Issue_of_Data
    read(10, '(3X, 4D19.12)') Transmission_time_of_message, Fit_interval, &
      spare1, spare2
    close(10)

    write(6, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2, F5.1,3D19.12)') &
      PRN, year, month, day, hour, minute, second, &
      clock_bias, clock_drift, clock_drift_rate
    write(6, '(3X, 4D19.12)') IODE, Crs, Delta_n, M0
    write(6, '(3X, 4D19.12)') Cuc, e, Cus, sqrt_A
    write(6, '(3X, 4D19.12)') TOE, Cic, LOMEGA0, Cis
    write(6, '(3X, 4D19.12)') i0, Crc, somega, OMEGA_DOT
    write(6, '(3X, 4D19.12)') IDOT, Codes_on_L2_channel, GPS_Week, &
      L2_P_data_flag
    write(6, '(3X, 4D19.12)') accuracy, health, TGD, IODC_Issue_of_Data
    write(6, '(3X, 4D19.12)') Transmission_time_of_message, Fit_interval, &
      spare1, spare2

  end subroutine read_GPS_Nav
end module file_io
