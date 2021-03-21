module navigation_message
  use time_util
  use exec_conditions
  implicit none

  ! 電離層遅延補正情報
  INTEGER, PARAMETER :: IONO_PARAMETERS = 4
  DOUBLE PRECISION :: ion_alpha(IONO_PARAMETERS)
  DOUBLE PRECISION :: ion_beta(IONO_PARAMETERS)

  ! うるう秒の情報
  INTEGER          :: leap_sec = 0

  ! エフェメリス情報を格納する構造体
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
  end type ephemeris_info

  TYPE(ephemeris_info) :: ephem_buf(MAX_PRN, MAX_EPHMS) ! 全エフェメリス格納配列
  TYPE(ephemeris_info) :: ephem_data ! 1衛星分のエフェメリス
  TYPE(ephemeris_info) :: current_ephem(MAX_PRN) ! 測位計算に使用するエフェメリス格納配列
  INTEGER              :: ephem_count(MAX_PRN) = 0

contains
  subroutine read_nav_msg()
    implicit none

    ! 使用局所領域
    INTEGER          :: ios ! ファイル読み込みステータス
    INTEGER          :: year, month, day, hour, minute ! TOC計算用
    DOUBLE PRECISION :: second ! TOC計算用
    INTEGER          :: i, j ! ループ用カウンタ
    DOUBLE PRECISION :: t ! (TOCの週番号 - GPSweek) + TOT [sec]
    type(wtime)      :: wt_toc
    INTEGER          :: info_week


    ! Navigation Message Fileオープン
    open(10, file=nav_msg_file, iostat=ios,action='read', status='old')
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

    do ! ファイルの最終行まで読み込む

      ! ----------- 1行目 読み込み ----------------------------------
      read(10, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12)', iostat = ios) &
        ephem_data%PRN, year, month, day, hour, minute, second, &
        ephem_data%AF0, ephem_data%AF1, ephem_data%AF2

      if (ios < 0) exit ! ファイル最終行に来たら読み込み終了

      wt_toc%week = 0 ! 週番号を初期化
      wt_toc%sec = 0.d0 ! 経過秒を初期化
      call date_to_wtime(year, month, day, hour, minute, second, wt_toc) ! TOCを計算する
      ephem_data%TOC = wt_toc%sec
      info_week = wt_toc%week


      ! ----------- 2行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%IODE, ephem_data%Crs, ephem_data%delta_n, ephem_data%M0
      ! ----------- 3行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%Cuc, ephem_data%e,  ephem_data%Cus, ephem_data%sqrtA
      ! ----------- 4行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%TOE, ephem_data%Cic, ephem_data%LOMEGA0, ephem_data%Cis
      ! ----------- 5行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%i0, ephem_data%Crc, ephem_data%somega, ephem_data%OMEGA_DOT
      ! ----------- 6行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%IDOT, ephem_data%CAonL2, ephem_data%WEEK, ephem_data%L2P
      ! ----------- 7行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)') ephem_data%acc, ephem_data%health, ephem_data%TGD, ephem_data%IODC
      ! ----------- 8行目 読み込み ----------------------------------
      read(10, '(3X, 4D19.12)', iostat=ios) ephem_data%TOT, ephem_data%Fit

      if (ephem_data%PRN < 1) cycle

      ! 週番号を揃える
      t = (info_week - ephem_data%WEEK) * WEEK_SEC
      ephem_data%TOC = ephem_data%TOC + t

      ! PRNが同じ衛星がまだひとつもないとき、エフェメリス格納配列に入れて、カウンタを1にする
      if (ephem_count(ephem_data%PRN) == 0) then
        ephem_buf(ephem_data%PRN, 1) = ephem_data ! i番目にデータを挿入
        ephem_count(ephem_data%PRN) = 1
        cycle ! 次のエフェメリスの読み込みにジャンプ
      end if

      ! 同じPRNのエフェメリスは送信時刻の早いものを残す-------------------------
      do  i = 1, ephem_count(ephem_data%PRN)
        if (ephem_buf(ephem_data%PRN, i)%WEEK /= wt_toc%week) cycle ! 週番号が異なるものは飛ばす
        if (ephem_buf(ephem_data%PRN, i)%WEEK == ephem_data%IODC) then !IODCが一致
          if (ephem_data%TOT < ephem_buf(ephem_data%PRN, i)%TOT) then ! 送信時刻が早いものを残す
            ephem_buf(ephem_data%PRN, i) = ephem_data
          end if
          ephem_data%PRN = 0
          exit
        end if
      end do
      ! ------------------------------------------------------------------

      if (ephem_data%PRN < 1) cycle

      if (ephem_count(ephem_data%PRN) >= MAX_EPHMS) & ! 同じPRNの衛星が多すぎるとエラー
          stop 'too long NAV file.'

      ! 配列に格納する(送信時刻の昇順) ----------------------------------------
      do i = 1, ephem_count(ephem_data%PRN)
        t = (ephem_buf(ephem_data%PRN, i)%WEEK - info_week) * WEEK_SEC &
            + ephem_buf(ephem_data%PRN, i)%TOT
        if (ephem_data%TOT < t) exit
      end do
      do j = ephem_count(ephem_data%PRN), i + 1, -1
        ! 現在のPRNの衛星のエフェメリスを送信時刻の昇順で上にずらし、i番目を開ける
        ephem_buf(ephem_data%PRN, j) = ephem_buf(ephem_data%PRN, j - 1)
      end do
      ephem_buf(ephem_data%PRN, i) = ephem_data ! i番目にデータを挿入
      ephem_count(ephem_data%PRN) = ephem_count(ephem_data%PRN) + 1 ! カウントを1増やす
      ! ------------------------------------------------------------------

    end do
    close(10)

  end subroutine read_nav_msg

! ============================================================================================

  subroutine set_ephemeris(PRN, wt, iode, return_code)
    implicit none
    ! 使用するエフェメリスをセットする
    ! 引数詳細
    INTEGER, INTENT(IN)                 :: PRN
    TYPE(wtime), INTENT(IN)             :: wt
    INTEGER, INTENT(IN)                 :: iode
    INTEGER, INTENT(INOUT)              :: return_code

    ! 使用局所領域
    DOUBLE PRECISION :: t0, t
    INTEGER          :: i  ! カウンタ

    ! リターンコードを初期化
    return_code = 0
    ! 新しいエフェメリスから探す
    do i = ephem_count(PRN), 1, -1
      ! j = i
      t0 = ( ephem_buf(PRN, i)%WEEK - wt%week ) * WEEK_SEC ! エフェメリスの週番号と指定時刻の週番号の差 [sec]
      t = ephem_buf(PRN, i)%TOC + t0 ! t0をクロックに足す
      if (wt%sec < t - EPHEMERIS_EXPIRE * 3600.d0 - 0.1d0 .or. & ! 有効期限より小さい
          wt%sec > t - EPHEMERIS_EXPIRE * 3600.d0 + 0.1d0) then  ! 有効期限より大きい
        cycle ! 有効期限内でなければ飛ばす
      end if

      ! IODEが一致すれば採用
      if (iode >= 0) then
        if (ephem_buf(PRN, i)%IODE == iode) exit
        cycle
      end if

      ! 受信時刻が指定時刻よりも前なら採用
      t = ephem_buf(PRN, i)%TOT + t0 ! tを送信時刻(TOT)に足す
      if (t < wt%sec + 0.1d0 ) exit

    end do
    if (i >= 0) then
      current_ephem(PRN) = ephem_buf(PRN, i)
    else
      return_code = 9
    end if
  end subroutine set_ephemeris


end module navigation_message
