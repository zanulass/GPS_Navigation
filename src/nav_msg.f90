module navigation_message
  use mod_variable
  use time_util
  implicit none
contains
  subroutine read_nav_msg(nav_msg_file)
    ! GPS NAVIGATION MESSAGE FILE からパラメータを読み込む

    CHARACTER(256), INTENT(IN)   :: nav_msg_file         ! RINEX NAVIGATION MESSAGE FILEのパス

    INTEGER        :: ios ! ファイル読み込みステータス
    INTEGER :: year, month, day, hour, minute ! TOC計算用
    DOUBLE PRECISION :: second ! TOC計算用
    INTEGER          :: i ! ループ用カウンタ
    type(wtime) :: wt

    ! ファイルオープン

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

      wt%week = 0 ! 週番号を初期化
      wt%sec = 0.d0 ! 経過秒を初期化
      call date_to_wtime(year, month, day, hour, minute, second, wt) ! TOCを計算する
      ephem_buf%TOC = wt%sec
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

  end subroutine read_nav_msg
end module navigation_message
