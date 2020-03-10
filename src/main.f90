program DO3SE_main

    use DO3SE_Photosynthesis_ml
    ! use mocked_functions
    use DO3SE_ConfigTypes_ml

    implicit none

    character(len=256) :: arg1

    !== Photosynthesis vars ==!
    type(pngstoconfig_t) :: pgc   !< Photosynthesis gsto parameters
    type(season_t) :: season
    real :: V_cmax_25   !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
    real :: J_max_25    !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)
    real :: D_0         !< "The VPD at which g_sto is reduced by a factor of 2" (kPa) (Leuning et al. 1998)
    real :: LAI
    real :: PARsun
    real :: PARshade
    real :: sinB
    real :: Lm          !< Leaf dimension (m)
    real :: Tair_C      !< Ambient air temperature (degrees C)
    real :: u           !< Wind speed (m/s)
    real :: CO2         !< CO2 concentration (ppm)
    real :: Rn          !< Net radiation (MJ m-2 h-1)
    real :: R           !< Global radiation (W m-2)
    real :: albedo      !< Surface albedo (fraction)
    real :: PPFD        !< PPFD (umol m-2 s-1)
    real :: P           !< Ambient air pressure (kPa)
    real :: eact        !< Ambient vapour pressure (kPa)
    real :: O3
    real :: td
    integer :: dd
    integer :: hr

    real, dimension(366) :: fO3_h_1_hist
    real, dimension(366) :: fO3_d_hist
    real, dimension(366) :: daily_thermal_temperatures
    real :: Tleaf_C  !< Leaf temperature (degrees C)
    real :: A_n        !< Output: Net CO2 assimilation (umol m-2 PLA s-1)
    real :: Canopy_A_n
    real :: A_c
    real :: A_j
    real :: A_p
    real :: O3up
    real :: O3up_acc
    real :: fO3_h
    real :: fO3_d
    real :: R_d
    real :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
    real :: g_sv       !< Output: Stomatal conductance to water vapour
    real :: g_bv       !< Output: Boundary conductance to water vapour
    !== END Photosynthesis vars ==!
  


    print *, "hello world"
    call get_command_argument(1, arg1)
    print * , arg1
    ! ============== TEST gsto_pn ============== !

    !== specify value from spanish wheat photo config ==!
    pgc%g_sto_0 = 20000.0
    pgc%m = 8.12
    pgc%v_cmax_25 = 180.0
    pgc%j_max_25 = 400.0
    
    season%growing_season_method = "constant"
    season%SGS = 118
    season%EGS = 210
    
    season%accumulation_period_method = "constant"
    season%Astart = 153
    season%Aend = 208
    
    ! season%height_method = "constant"
    ! season%height = 1.0
    
    ! season%LAI_method = "day PLF"
    season%LAI_a = 0.0
    season%LAI_b = 3.5
    season%LAI_c = 3.5
    season%LAI_d = 0.0
    season%LAI_1 = 21
    season%LAI_2 = 21
    
    season%SAI_method = "wheat"

    V_cmax_25 = 180.0   !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
    J_max_25 = 400.0    !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)
    D_0 = 1000       !< "The VPD at which g_sto is reduced by a factor of 2" (kPa) (Leuning et al. 1998)
    LAI = 0
    PARsun = 1
    PARshade = 1
    sinB = 1
    Lm = 0.02         !< Leaf dimension (m)
    Tair_C = 20      !< Ambient air temperature (degrees C)
    u = 1           !< Wind speed (m/s)
    CO2 = 1         !< CO2 concentration (ppm)
    Rn = 1          !< Net radiation (MJ m-2 h-1)
    R = 1           !< Global radiation (W m-2)
    albedo = 0.2     !< Surface albedo (fraction)
    PPFD = 1        !< PPFD (umol m-2 s-1)
    P = 1           !< Ambient air pressure (kPa)
    eact = 1        !< Ambient vapour pressure (kPa)
    O3 = 50
    td = 1
    dd = 100
    hr = 2
  
      
    call gsto_pn(pgc, season, V_cmax_25, J_max_25, D_0, &
    LAI, PARsun, PARshade, sinB, &
    Lm, Tair_C, u, CO2, &
    PPFD, Rn, R, albedo, P, eact, O3, td, dd, hr, &
    fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
    Tleaf_C, A_n, Canopy_A_n, A_c, A_j, A_p, O3up, &
    O3up_acc, fO3_h, fO3_d, R_d, g_sto, g_sv, g_bv)

    ! call gsto_pn(g_sto)
    print *, g_sto
  
end program DO3SE_main
  