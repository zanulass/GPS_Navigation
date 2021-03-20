module correction_parameters
    use mod_variable
contains
  subroutine correct_sat_clock(wt, current_ephem, sat_clock)
    ! 衛星のクロック誤差を修正
    ! 引数詳細
    TYPE(wtime), INTENT(IN) :: wt
    TYPE(ephemeris_info), INTENT(IN) :: current_ephem
    DOUBLE PRECISION, INTENT(OUT) :: sat_clock

   ! 使用局所領域
    DOUBLE PRECISION :: tk0 = 0.d0, tk = 0.d0
    DOUBLE PRECISION :: dt = 0.d0! クロックの 補正量
    DOUBLE PRECISION :: tr = 0.d0

    ! 時間差を求める
    tk0 = ( wt%week - current_ephem%WEEK ) * WEEK_SEC + wt%sec - current_ephem%TOC
    tk = tk0

    ! 衛星時計の補正量を計算
    dt = current_ephem%AF0 &
       + current_ephem%AF1 * tk &
       + current_ephem%AF2 * tk * tk

    sat_clock = dt + tr - current_ephem%TGD
  end subroutine correct_sat_clock

end module

