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
    INTEGER          :: i, j ! ループ用カウンタ
    DOUBLE PRECISION :: t ! (TOCの週番号 - GPSweek) + TOT [sec]
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

    do ! ファイルの最終行まで読み込む

      ! ----------- 1行目 読み込み ----------------------------------
      read(10, '(I2,1X,I2.2,1X,I2,1X,I2,1X,I2,1X,I2,F5.1,3D19.12)', iostat = ios) &
        ephem_data%PRN, year, month, day, hour, minute, second, &
        ephem_data%AF0, ephem_data%AF1, ephem_data%AF2

      if (ios < 0) exit ! ファイル最終行に来たら読み込み終了

      wt%week = 0 ! 週番号を初期化
      wt%sec = 0.d0 ! 経過秒を初期化
      call date_to_wtime(year, month, day, hour, minute, second, wt) ! TOCを計算する
      ephem_data%TOC = wt%sec
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

      ! 週番号を揃える
      t = (wt%week - ephem_data%WEEK) * WEEK_SEC
      ephem_data%TOC = ephem_data%TOC + t

      ! PRNが同じ衛星がまだひとつもないとき、エフェメリス格納配列に入れて、カウンタを1増やす
      if (ephem_count(ephem_data%PRN) == 0) then
        ephem_buf(ephem_data%PRN, 1) = ephem_data ! i番目にデータを挿入
        ephem_count(ephem_data%PRN) = 1
        cycle ! 次のエフェメリスの読み込みにジャンプ
      end if

      ! 同じPRNのエフェメリスは送信時刻の早いものを残す-------------------------
      ! ephem_bufにインデックス0が無いため、すでにephem_bufに1つ以上データが入っている場合に分岐に入る
      if (ephem_count(ephem_data%PRN) > 0) then
        do  i = 1, ephem_count(ephem_data%PRN)
          if (ephem_buf(ephem_data%PRN, i)%WEEK /= wt%week) cycle ! 週番号が同じ
          if (ephem_buf(ephem_data%PRN, i)%WEEK /= ephem_data%IODC) cycle !IODCが同じ
          if (ephem_data%TOT < ephem_buf(ephem_data%PRN, i)%TOT) then ! 送信時刻が早いものを残す
            ephem_buf(ephem_data%PRN, i) = ephem_data
          end if
        end do
      end if
      ! ------------------------------------------------------------------


      if (ephem_count(ephem_data%PRN) >= MAX_EPHMS) & ! 同じPRNの衛星が多すぎるとエラー
          stop 'too long NAV file.'

      ! 配列に格納する(送信時刻の昇順) ----------------------------------------
      do i = 1, ephem_count(ephem_data%PRN) + 1
          t = (ephem_buf(ephem_data%PRN, i)%WEEK - ephem_data%WEEK) * WEEK_SEC &
              + ephem_buf(ephem_data%PRN, i)%TOT
          ! 送信時刻がi番目のエフェメリスより早いとき、i番目に現在のエフェメリスを挿入する
          if (ephem_data%TOT < t) then
            do j = ephem_count(ephem_data%PRN) + 1, i, -1
              ! 現在のPRNの衛星のエフェメリスを送信時刻の昇順で上にずらし、i番目を開ける
              ephem_buf(ephem_data%PRN, j + 1) = ephem_buf(ephem_data%PRN, j)
            end do
            ephem_buf(ephem_data%PRN, i) = ephem_data ! i番目にデータを挿入
            exit
          end if
      end do
      ephem_count(ephem_data%PRN) = ephem_count(ephem_data%PRN) + 1 ! カウントを1増やす
      ! ------------------------------------------------------------------

    end do
    close(10)

  end subroutine read_nav_msg

  subroutine set_ephemeris(PRN, wt, iode, current_ephem)
    ! 使用するエフェメリスをセットする
    ! 引数詳細
    INTEGER, INTENT(IN) :: PRN
    TYPE(wtime), INTENT(IN) :: wt
    INTEGER, INTENT(IN) :: iode
    TYPE(ephemeris_info), INTENT(INOUT) :: current_ephem

    ! 使用局所領域
    DOUBLE PRECISION :: t0, t
    INTEGER :: i, j ! カウンタ

    j = 0
    ! 新しいエフェメリスから探す
    do i = ephem_count(PRN), 1, -1
      j = i
      ! エフェメリスの週番号と指定時刻の週番号の差 [sec]
      t0 = ( ephem_buf(PRN, i)%WEEK - wt%week ) * WEEK_SEC
      t = ephem_buf(PRN, i)%TOC + t0 ! tをクロックに足す
      if (wt%sec < t - EPHEMERIS_EXPIRE * 3600.d0 - 0.1d0 .or. & ! 有効期限より小さい
          wt%sec > t - EPHEMERIS_EXPIRE * 3600.d0 + 0.1d0) then  ! 有効期限より大きい
        cycle ! 有効期限内でなければ飛ばす
      end if

      ! IODEが一致すれば採用
      if (iode >= 0) then
        if (ephem_buf(PRN, i)%IODE == iode) then
          j = i
          exit
        end if
      end if

      ! 受信時刻が指定時刻よりも前なら採用
      t = ephem_buf(PRN, i)%TOT + t0 ! tを送信時刻(TOT)に足す
      if (t < wt%sec + 0.1d0 ) then
        j = i
        exit
      end if

    end do

    current_ephem = ephem_buf(PRN, j)
  end subroutine set_ephemeris


end module navigation_message
