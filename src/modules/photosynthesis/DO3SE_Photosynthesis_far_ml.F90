module DO3SE_Photosynthesis_far_ml
    ! #define TESTMODE
    ! #ifdef TESTMODE
    !     use DO3SE_Photosynthesis_mocks
    ! #else
        use DO3SE_Photosynthesis_mocks
        use DO3SE_Photosynthesis_helpers_ml
        use DO3SE_PhysicalConstants_ml, only: T0, R, Rn_MJ_to_W, DRATIO
        use DO3SE_ModelConstants_ml, only: UNDEF, LEAF_G_H2O
        ! use DO3SE_Met_ml, only: sunlit_LAI
        use DO3SE_ConfigTypes_ml
        use DO3SE_Util_ml
    ! #endif

    implicit none
    private

    ! public :: farquhar_photosynthesis
    public :: farquhar_photosynthesis_2

    contains

        !> Model the net assimilation rate of C3 plants according to the model
        !! developed by Farquar (1980).  Outputs are stomatal conductance and net
        !! CO2 assimilation.
        ! pure subroutine farquhar_photosynthesis(Tleaf_C, c_a, e_a, Q, g_bl, &
        !                                         g_sto_0, m, V_cmax_25, J_max_25, D_0, &
        !                                         gsto_final, pngsto_An, pngsto_Ac, pngsto_Aj, pngsto_Ap, pngsto_R_d)
        ! real, intent(in) :: Tleaf_C       !< Leaf temperature (degrees C)
        ! real, intent(in) :: c_a           !< CO2 concentration (ppm)
        ! real, intent(in) :: e_a           !< Ambient vapour pressure (Pa)
        ! real, intent(in) :: Q             !< PPFD (umol/m^2/s)
        ! real, intent(in) :: g_bl          !< Boundary layer conductance to H2O vapour (micromol m-2 PLA s-1)
        ! real, intent(in) :: g_sto_0       !< Closed stomata conductance (umol/m^2/s)
        ! real, intent(in) :: m             !< Species-specific sensitivity to An (dimensionless)
        ! real, intent(in) :: V_cmax_25     !< Maximum catalytic rate at 25 degrees (umol/m^2/s)
        ! real, intent(in) :: J_max_25      !< Maximum rate of electron transport at 25 degrees (umol/m^2/s)
        ! real, intent(in) :: D_0           !< "The VPD at which g_sto is reduced by a factor of 2" (Pa) (Leuning et al. 1998)

        ! real, intent(out) :: gsto_final   !< Output: Raw photosynthesis-based stomatal conductance (umol m-2 s-1)
        ! real, intent(out) :: pngsto_An    !< Output: net CO2 assimilation (umol m-2 s-1)
        ! real, intent(out) :: pngsto_Ac
        ! real, intent(out) :: pngsto_Aj
        ! real, intent(out) :: pngsto_Ap
        ! real, intent(out) :: pngsto_R_d

        ! ! parameters considered (or defined) to be constant for all species
        ! real, parameter :: O_i = 210.0           !O2 concentration                   [mmol/mol]
        ! real, parameter :: E_K_C = 79430.0       !activation energy of K_C           [J/mol]            Medlyn2002
        ! real, parameter :: E_K_O = 36380.0       !activation energy of K_O           [J/mol]            Medlyn2002
        ! real, parameter :: E_R_d = 53000.0       !activation energy of R_d           [J/mol]            Leuning1995
        ! real, parameter :: E_Gamma_star = 37830.0 !activation energy for C-comp-point [J/mol]            Medlyn2002
        ! real, parameter :: K_C_25 = 404.9        !K.C at reference temperature 25    [micro mol/mol]    Medlyn2002
        ! real, parameter :: K_O_25 = 278.4        !K.O at reference temperature 25    [mmol/mol]         Medlyn2002
        ! real, parameter :: R_d_20 = 0.32         !R_d at reference temperature 20    [micro mol/(m^2*s)]Leuning1995
        ! real, parameter :: Gamma_star_25 = 42.75 !CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002
        ! real, parameter :: A_j_a = 4.0           !electron requirement for NADPH formation
        ! real, parameter :: A_j_b = 8.0           !electron requirement for ATP formation

        ! ! species spedific model parameters (that don't tend to have species specific
        ! ! values, others are supplied as arguments)
        ! real, parameter :: alpha = 0.3           !efficiency light energy conversion [mol electrons/mol photons]
        ! real, parameter :: Teta = 0.90           !shape of J~Q determining factor    []
        ! real, parameter :: H_a_jmax = 50300      !activation energy for J_max        [J/mol]
        ! real, parameter :: H_d_jmax = 152044     !deactivation energy for J_max      [J/mol]
        ! real, parameter :: H_a_vcmax = 73637     !activation energy for V_cmax       [J/mol]
        ! real, parameter :: H_d_vcmax = 149252    !deactivation energy for V_cmax     [J/mol]
        ! real, parameter :: S_V_vcmax = 486       !entropy terms                      [J/(mol*K)]
        ! real, parameter :: S_V_jmax = 495        !entropy terms                      [J/(mol*K)

        ! ! Converted inputs
        ! real :: Tleaf_K

        ! ! state variables
        ! real :: A_n                             !netto assimilation rate            [micro mol/(m^2*s)]
        ! real :: A_c                             !Rub. activity. lim. ass. rate      [micro mol/(m^2*s)]
        ! real :: A_j                             !electr. transp. lim. ass. rate     [micro mol/(m^2*s)]
        ! real :: A_p                             !triose phosphate utilisation lim. ass. rate   [micro mol/(m^2*s)]
        ! real :: Gamma_star                      !CO2 comp. point without day resp.  [micro mol/mol]
        ! real :: R_d                             !day respiration rate               [micro mol/(m^2*s)]
        ! real :: K_C                             !Michaelis constant CO2             [micro mol/mol]
        ! real :: K_O                             !Michaelis constant O2              [mmol/mol]
        ! real :: J                               !Rate of electron transport         [micro mol/(m^2*s)]
        ! real :: Gamma                           !CO2 compensation point             [micro mol/mol]
        ! real :: J_max                           !Max rate of electron transport     [micro mol/(m^2*s)]
        ! real :: V_cmax                          !Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
        ! real :: e_sat_i                         !internal saturation vapour pressure[Pa]
        ! real :: g_sto                           !two sided stomatal conduct.,vapour [micro mol/(m^2s)]
        ! real :: h_s                             !relative humidity at leaf surface  [decimal fraction]
        ! real :: h_s_VPD                         !VPD at leaf surface                [Pa]
        ! real :: c_s                             !CO2 concentration at leaf surface  [micromol/mol]
        ! real :: c_i                             !CO2 concentration inside stomata   [micromol/mol]

        ! ! iteration parameters
        ! integer :: iterations                   !number of the iterations bofore convergence
        ! real :: c_i_sup                         !CO2 concentration inside stomata possible through supply
        ! integer :: k                            !loop parameters

        ! Tleaf_K = T0 + Tleaf_C

        ! ! Calculation of the model variables which are only
        ! ! dependend on environmental conditions:

        ! Gamma_star = temp_dep(Gamma_star_25, T0 + 25, E_Gamma_star, Tleaf_K, R)

        ! K_C = temp_dep(K_C_25, T0 + 25, E_K_C, Tleaf_K, R)

        ! K_O = temp_dep(K_O_25, T0 + 25, E_K_O, Tleaf_K, R)

        ! R_d = temp_dep(R_d_20, T0 + 20, E_R_d, Tleaf_K, R)

        ! J_max = temp_dep_inhibit(J_max_25, T0 + 25, H_a_jmax, H_d_jmax, S_V_jmax, Tleaf_K, R)

        ! V_cmax = temp_dep_inhibit(V_cmax_25, T0 + 25, H_a_vcmax, H_d_vcmax, S_V_vcmax, Tleaf_K, R)

        ! ! Electron transport rate
        ! J = (J_max + alpha*Q - sqrt((J_max + alpha*Q)**2 - 4*alpha*Q*Teta*J_max)) / (2*Teta)

        ! e_sat_i    = 1000 * saturated_vapour_pressure(Tleaf_C)

        ! Gamma        = (Gamma_star+(K_C*R_d*(1+(O_i/K_O))/V_cmax))/&
        !                 (1-(R_d/V_cmax))

        ! !The following loop guesses a start value for c_i and tests whether
        ! !it satisfies all the relevant restrictions. If not a new value for
        ! !c_i is tested:

        ! c_i         = 0.0

        ! ! gsto needs a starting point, so let's set it to g_sto_0
        ! g_sto = g_sto_0

        ! do k=1,50

        !     ! Rubisco activity limited assimilation rate
        !     A_c = V_cmax * ((c_i - Gamma_star) &
        !                     / (c_i + (K_C * (1 + (O_i / K_O)))))

        !     ! RuBP regeneration (electron transport) limited assimilation rate
        !     A_j = J * ((c_i - Gamma_star) &
        !             / ((A_j_a * c_i) + (A_j_b * Gamma_star)))

        !     ! Triose phosphate utilisation limited assimilation rate
        !     A_p = 0.5 * V_cmax

        !     ! CO2 assimilation rate
        !     A_n = min(A_c, A_j, A_p) - R_d

        !     ! Surface CO2
        !     c_s = c_a - (A_n * (1.37/g_bl))

        !     ! Surface humidity
        !     h_s = (g_sto*e_sat_i + g_bl*e_a) / (e_sat_i * (g_sto + g_bl))
        !     ! Convert relative humidity to VPD
        !     h_s_VPD = e_sat_i - (e_sat_i * h_s)

        !     ! Stomatal conductance
        !     ! TODO: use humidity deficit version instead
        !     g_sto = g_sto_0 + m * (A_n / ((1 + (h_s_VPD/D_0)) * (c_s - Gamma)))*1e6

        !     ! CO2 supply
        !     c_i_sup = c_a - ((A_n*(1.6/g_sto + 1.37/g_bl))*1e6)

        !     !exits the loop when c_i calculated with both ways meet the convergence
        !     !criterium:

        !     iterations = k
        !     if (abs(c_i - c_i_sup) < 0.001) then
        !     exit
        !     end if

        !     !Guesses a new c_i as the mean of the first guess and c_i resulting from
        !     !the supply function:

        !     c_i      = c_i-(c_i-c_i_sup)/2

        ! end do

        ! ! Calculate final stomatal conductances
        ! gsto_final = g_sto
        ! pngsto_An = A_n
        ! pngsto_Ac = A_c
        ! pngsto_Aj = A_j
        ! pngsto_Ap = A_p
        ! pngsto_R_d = R_d
        ! end subroutine farquhar_photosynthesis


        pure subroutine farquhar_photosynthesis_2(pgc, season, Tleaf_C, c_a, e_a, Q, g_bl, &
                                                    g_sto_0, m, V_cmax_25, J_max_25, D_0, O3, td, dd, hr, &
                                                    fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
                                                    gsto_final, pngsto_An, pngsto_Ac, &
                                                    pngsto_Aj, pngsto_Ap, O3up_out, O3up_acc_out, &
                                                    fO3_h_out, fO3_d_out, pngsto_R_d)

            type(pngstoconfig_t), intent(in) :: pgc   !< Photosynthesis gsto parameters
            type(season_t), intent(in) :: season
            real, intent(in) :: Tleaf_C       !< Leaf temperature (degrees C)
            real, intent(in) :: c_a           !< CO2 concentration (ppm)
            real, intent(in) :: e_a           !< Ambient vapour pressure (Pa)
            real, intent(in) :: Q             !< PPFD (umol/m^2/s)
            real, intent(in) :: g_bl          !< Boundary layer conductance to H2O vapour (micromol m-2 PLA s-1)
            real, intent(in) :: g_sto_0       !< Closed stomata conductance (umol/m^2/s)
            real, intent(in) :: m             !< Species-specific sensitivity to An (dimensionless)
            real, intent(in) :: V_cmax_25     !< Maximum catalytic rate at 25 degrees (umol/m^2/s)
            real, intent(in) :: J_max_25      !< Maximum rate of electron transport at 25 degrees (umol/m^2/s)
            real, intent(in) :: D_0           !< "The VPD at which g_sto is reduced by a factor of 2" (Pa) (Leuning et al. 1998)
            real, intent(in) :: O3
            real, intent(in) :: td
            integer, intent(in) :: dd
            integer, intent(in) :: hr

            real, dimension(366), intent(inout) :: fO3_h_1_hist
            real, dimension(366), intent(inout) :: fO3_d_hist
            real, dimension(366), intent(inout) :: daily_thermal_temperatures
            real, intent(out) :: gsto_final   !< Output: Raw photosynthesis-based stomatal conductance (umol m-2 s-1)
            real, intent(out) :: pngsto_An    !< Output: net CO2 assimilation (umol m-2 s-1)
            real, intent(out) :: pngsto_Ac
            real, intent(out) :: pngsto_Aj
            real, intent(out) :: pngsto_Ap
            real, intent(out) :: pngsto_R_d
            real, intent(out) :: O3up_out
            real, intent(inout) :: O3up_acc_out
            real, intent(out) :: fO3_h_out
            real, intent(out) :: fO3_d_out

            ! parameters considered (or defined) to be constant for all species
            real, parameter :: O_i = 210.0           !O2 concentration                   [mmol/mol]
            real, parameter :: E_K_C = 79430.0       !activation energy of K_C           [J/mol]            Medlyn2002
            real, parameter :: E_K_O = 36380.0       !activation energy of K_O           [J/mol]            Medlyn2002
            real, parameter :: E_R_d = 53000.0       !activation energy of R_d           [J/mol]            Leuning1995
            real, parameter :: E_Gamma_star = 37830.0 !activation energy for C-comp-point [J/mol]            Medlyn2002
            real, parameter :: K_C_25 = 404.9        !K.C at reference temperature 25    [micro mol/mol]    Medlyn2002
            real, parameter :: K_O_25 = 278.4        !K.O at reference temperature 25    [mmol/mol]         Medlyn2002
            real, parameter :: R_d_20 = 0.32         !R_d at reference temperature 20    [micro mol/(m^2*s)]Leuning1995
            real, parameter :: Gamma_star_25 = 42.75 !CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002
            real, parameter :: A_j_a = 4.0           !electron requirement for NADPH formation
            real, parameter :: A_j_b = 8.0           !electron requirement for ATP formation

            ! species spedific model parameters (that don't tend to have species specific
            ! values, others are supplied as arguments)
            real, parameter :: alpha = 0.3           !efficiency light energy conversion [mol electrons/mol photons]
            real, parameter :: Teta = 0.90           !shape of J~Q determining factor    []
            real, parameter :: H_a_jmax = 50300      !activation energy for J_max        [J/mol]
            real, parameter :: H_d_jmax = 152044     !deactivation energy for J_max      [J/mol]
            real, parameter :: H_a_vcmax = 73637     !activation energy for V_cmax       [J/mol]
            real, parameter :: H_d_vcmax = 149252    !deactivation energy for V_cmax     [J/mol]
            real, parameter :: S_V_vcmax = 486       !entropy terms                      [J/(mol*K)]
            real, parameter :: S_V_jmax = 495        !entropy terms                      [J/(mol*K)
            real, parameter :: fDO3 = 0.93
            real, parameter :: Gamma_3 = 0.5

            ! Converted inputs
            real :: Tleaf_K

            ! state variables
            real :: A_n                             !netto assimilation rate            [micro mol/(m^2*s)]
            real :: A_c                             !Rub. activity. lim. ass. rate      [micro mol/(m^2*s)]
            real :: A_j                             !electr. transp. lim. ass. rate     [micro mol/(m^2*s)]
            real :: A_p                             !triose phosphate utilisation lim. ass. rate   [micro mol/(m^2*s)]
            real :: O3up
            real :: O3up_acc
            real :: fO3_h
            real :: fO3_d
            real :: fO3_l
            real :: f_LA
            real :: f_LS
            real :: Gamma_star                      !CO2 comp. point without day resp.  [micro mol/mol]
            real :: R_d                             !day respiration rate               [micro mol/(m^2*s)]
            real :: K_C                             !Michaelis constant CO2             [micro mol/mol]
            real :: K_O                             !Michaelis constant O2              [mmol/mol]
            real :: J                               !Rate of electron transport         [micro mol/(m^2*s)]
            real :: Gamma                           !CO2 compensation point             [micro mol/mol]
            real :: J_max                           !Max rate of electron transport     [micro mol/(m^2*s)]
            real :: V_cmax                          !Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
            real :: e_sat_i                         !internal saturation vapour pressure[Pa]
            real :: g_sto                           !two sided stomatal conduct.,vapour [micro mol/(m^2s)]
            real :: h_s                             !relative humidity at leaf surface  [decimal fraction]
            real :: h_s_VPD                         !VPD at leaf surface                [Pa]
            real :: c_s                             !CO2 concentration at leaf surface  [micromol/mol]
            real :: c_i                             !CO2 concentration inside stomata   [micromol/mol]
            real :: t_l
            real :: t_lem
            real :: t_lma
            real :: t_lep
            real :: t_lse
            real :: rO3_s
            real :: td_dd

            ! iteration parameters
            integer :: iterations                   !number of the iterations bofore convergence
            real :: c_i_sup                         !CO2 concentration inside stomata possible through supply
            integer :: k                            !loop parameters

            if (hr == 0) then
            daily_thermal_temperatures(dd) = td
            end if

            if (dd >= season%Astart) then
            td_dd = td - daily_thermal_temperatures(season%Astart)
            else
            td_dd = 0
            end if

            t_l = pgc%t_l
            t_lem = pgc%t_lem
            t_lma = pgc%t_lma
            t_lep = pgc%t_lep
            t_lse = pgc%t_lse

            Tleaf_K = T0 + Tleaf_C

            ! Calculation of the model variables which are only
            ! dependend on environmental conditions:

            Gamma_star = temp_dep(Gamma_star_25, T0 + 25, E_Gamma_star, Tleaf_K, R)

            K_C = temp_dep(K_C_25, T0 + 25, E_K_C, Tleaf_K, R)

            K_O = temp_dep(K_O_25, T0 + 25, E_K_O, Tleaf_K, R)

            R_d = temp_dep(R_d_20, T0 + 20, E_R_d, Tleaf_K, R)

            J_max = temp_dep_inhibit(J_max_25, T0 + 25, H_a_jmax, H_d_jmax, S_V_jmax, Tleaf_K, R)

            V_cmax = temp_dep_inhibit(V_cmax_25, T0 + 25, H_a_vcmax, H_d_vcmax, S_V_vcmax, Tleaf_K, R)

            ! Electron transport rate
            J = (J_max + alpha*Q - sqrt((J_max + alpha*Q)**2 - 4*alpha*Q*Teta*J_max)) / (2*Teta)

            e_sat_i    = 1000 * saturated_vapour_pressure(Tleaf_C)

            Gamma        = (Gamma_star+(K_C*R_d*(1+(O_i/K_O))/V_cmax))/&
            (1-(R_d/V_cmax))

            !The following loop guesses a start value for c_i and tests whether
            !it satisfies all the relevant restrictions. If not a new value for
            !c_i is tested:

            c_i         = 0.0

            ! gsto needs a starting point, so let's set it to g_sto_0
            g_sto = g_sto_0

            do k=1,50

            O3up = O3 * g_sto * fDO3

            fO3_h = calc_fO3_h(O3up)

            f_LA = calc_f_LA(t_lem, t_lma, t_l, td_dd)

            if (dd == 1) then
            rO3_s = f_LA
            else
            rO3_s = fO3_d_hist(dd-1) + (1 - fO3_d_hist(dd-1)) * f_LA
            end if

            if (dd == 1) then
            fO3_d = rO3_s
            else
            fO3_d = fO3_h_1_hist(dd) * rO3_s
            end if

            O3up_acc = O3up_acc_out + O3up

            fO3_l = 1 - (Gamma_3 * O3up_acc_out) ! USING YESTERDAY'S VALUE, COULD CHANGE TO ABOVE.

            t_lma = (t_lep + t_lse) * fO3_l

            t_lse = 0.33 * t_lma

            f_LS = calc_f_LS(t_lem, t_lep, t_lma, t_l, fO3_l, td_dd)

            ! Rubisco activity limited assimilation rate
            A_c = V_cmax * ((c_i - Gamma_star) &
            / (c_i + (K_C * (1 + (O_i / K_O))))) * fO3_d * f_LS

            ! RuBP regeneration (electron transport) limited assimilation rate
            A_j = J * ((c_i - Gamma_star) &
            / ((A_j_a * c_i) + (A_j_b * Gamma_star)))

            ! Triose phosphate utilisation limited assimilation rate
            A_p = 0.5 * V_cmax

            ! CO2 assimilation rate
            A_n = min(A_c, A_j, A_p) - R_d

            ! Surface CO2
            c_s = c_a - (A_n * (1.37/g_bl))

            ! Surface humidity
            h_s = (g_sto*e_sat_i + g_bl*e_a) / (e_sat_i * (g_sto + g_bl))
            ! Convert relative humidity to VPD
            h_s_VPD = e_sat_i - (e_sat_i * h_s)

            ! Stomatal conductance
            ! TODO: use humidity deficit version instead
            g_sto = g_sto_0 + m * (A_n / ((1 + (h_s_VPD/D_0)) * (c_s - Gamma)))*1e6

            ! CO2 supply
            c_i_sup = c_a - ((A_n*(1.6/g_sto + 1.37/g_bl))*1e6)

            !exits the loop when c_i calculated with both ways meet the convergence
            !criterium:

            iterations = k
            if (abs(c_i - c_i_sup) < 0.001) then
            exit
            end if

            !Guesses a new c_i as the mean of the first guess and c_i resulting from
            !the supply function:

            c_i      = c_i-(c_i-c_i_sup)/2

            end do

            ! Calculate final stomatal conductances
            fO3_d_hist(dd) = fO3_d
            if (hr == 0) then
            fO3_h_1_hist(dd) = fO3_h
            end if

            gsto_final = g_sto
            pngsto_An = A_n
            pngsto_Ac = A_c
            pngsto_Aj = A_j
            pngsto_Ap = A_p
            pngsto_R_d = R_d
            O3up_out = O3up
            O3up_acc_out = O3up_acc_out + O3up
            fO3_h_out = fO3_h
            fO3_d_out = fO3_d

        end subroutine farquhar_photosynthesis_2

end module DO3SE_Photosynthesis_far_ml