module DO3SE_Photosynthesis_ml

    use DO3SE_Photosynthesis_mocks
    use DO3SE_Photosynthesis_helpers_ml
    use DO3SE_Photosynthesis_far_ml
    use DO3SE_PhysicalConstants_ml, only: T0, R, Rn_MJ_to_W, DRATIO
    use DO3SE_ModelConstants_ml, only: UNDEF, LEAF_G_H2O
!     use DO3SE_Met_ml, only: sunlit_LAI
    use DO3SE_ConfigTypes_ml
    use DO3SE_Util_ml
!   #include "DO3SE_Util_ml.h"
    ! use DO3SE_Met_ml, only: saturated_vapour_pressure, &
    !                         leaf_net_radiation, &
    !                         leaf_temp_EB, &
    !                         leaf_temp_de_Boeck
    use DO3SE_Resistance_ml, only: leaf_gb

    implicit none
    private

    public :: gsto_pn
    ! public :: nikolov_leaf_temperature

  contains

    !> Use Farquar photosynthesis model to calculate: net CO2 assimilation,
    !! A_n (umol CO2 m-2 PLA s-1); stomatal conductance, g_sto (mmol O3 m-2 PLA s-1);
    !! and potentially leaf temperature, Tleaf_C (degrees C) if configured to do so.
    pure subroutine gsto_pn(pgc, season, V_cmax_25, J_max_25, D_0, &
                            LAI, PARsun, PARshade, sinB, &
                            Lm, Tair_C, u, CO2, &
                            PPFD, Rn, R, albedo, P, eact, O3, td, dd, hr, &
                            fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
                            Tleaf_C, A_n, Canopy_A_n, A_c, A_j, A_p, O3up, &
                            O3up_acc, fO3_h, fO3_d, R_d, g_sto, g_sv, g_bv)
      type(pngstoconfig_t), intent(in) :: pgc   !< Photosynthesis gsto parameters
      type(season_t), intent(in) :: season
      real, intent(in) :: V_cmax_25   !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
      real, intent(in) :: J_max_25    !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)
      real, intent(in) :: D_0         !< "The VPD at which g_sto is reduced by a factor of 2" (kPa) (Leuning et al. 1998)
      real, intent(in) :: LAI
      real, intent(in) :: PARsun
      real, intent(in) :: PARshade
      real, intent(in) :: sinB
      real, intent(in) :: Lm          !< Leaf dimension (m)
      real, intent(in) :: Tair_C      !< Ambient air temperature (degrees C)
      real, intent(in) :: u           !< Wind speed (m/s)
      real, intent(in) :: CO2         !< CO2 concentration (ppm)
      real, intent(in) :: Rn          !< Net radiation (MJ m-2 h-1)
      real, intent(in) :: R           !< Global radiation (W m-2)
      real, intent(in) :: albedo      !< Surface albedo (fraction)
      real, intent(in) :: PPFD        !< PPFD (umol m-2 s-1)
      real, intent(in) :: P           !< Ambient air pressure (kPa)
      real, intent(in) :: eact        !< Ambient vapour pressure (kPa)
      real, intent(in) :: O3
      real, intent(in) :: td
      integer, intent(in) :: dd
      integer, intent(in) :: hr

      real, dimension(366), intent(inout) :: fO3_h_1_hist
      real, dimension(366), intent(inout) :: fO3_d_hist
      real, dimension(366), intent(inout) :: daily_thermal_temperatures
      real, intent(inout) :: Tleaf_C  !< Leaf temperature (degrees C)
      real, intent(out) :: A_n        !< Output: Net CO2 assimilation (umol m-2 PLA s-1)
      real, intent(out) :: Canopy_A_n
      real, intent(out) :: A_c
      real, intent(out) :: A_j
      real, intent(out) :: A_p
      real, intent(out) :: O3up
      real, intent(out) :: O3up_acc
      real, intent(out) :: fO3_h
      real, intent(out) :: fO3_d
      real, intent(out) :: R_d
      real, intent(out) :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
      real, intent(out) :: g_sv       !< Output: Stomatal conductance to water vapour
      real, intent(out) :: g_bv       !< Output: Boundary conductance to water vapour

      integer, parameter :: MAX_ITERATIONS = 10
      integer, parameter :: MAX_TLEAF_ITERATIONS = 5

      integer :: i, j
      real :: R_ni
      real :: A_n_sun, A_n_shade, LAIsun, LAIshade

      ! aproximates the boundary layer conductance for forced convection
      ! (converted to umol m-2 s-1)
      g_bv = leaf_gb(LEAF_G_H2O, Lm, max(0.01, u)) * 1e6

      select case (pgc%Tleaf_method)
      ! case ("input")
      !   ! Tleaf_C supplied, just run photosynthesis
      !   call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                g_sv, A_n, A_c, A_j, A_p, R_d)
      ! case ("ambient")
      !   ! Use ambient air temperature as leaf temperature
      !   Tleaf_C = Tair_C
      !   call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                g_sv, A_n, A_c, A_j, A_p, R_d)

      ! case ("Nikolov")
      !   ! Estimate Tleaf_C according to Nikolov (1995), iteratively solved with
      !   ! A_n and g_sto, starting with Tleaf_C = Tair_C.
      !   Tleaf_C = Tair_C
      !   do i = 1, MAX_ITERATIONS
      !     call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                  pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                  g_sv, A_n, A_c, A_j, A_p, R_d)
      !     ! TODO: leaf wetness status
      !     Tleaf_C = nikolov_leaf_temperature(Tair_C, P*1e3, Rn*Rn_MJ_to_W, eact*1e3, &
      !                                        g_sv, g_bv, 0.0)
      !   end do
      ! case ("EB")
      !   ! Estimate Tleaf_C according to Campbell & Norman (1998), iteratively
      !   ! solved with A_n and g_sto, starting with Tleaf_C = Tair_C.
      !   Tleaf_C = Tair_C
      !   call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                g_sv, A_n, A_c, A_j, A_p, R_d)
      !   do i = 1, MAX_ITERATIONS
      !     do j = 1, MAX_TLEAF_ITERATIONS
      !       R_ni = leaf_net_radiation(R, eact, Tair_C, Tleaf_C, albedo, 1.0)
      !       Tleaf_C = leaf_temp_EB(Lm, Tair_C, P, eact, R_ni, u, g_sv*1e-6)
      !     end do
      !     call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                  pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                  g_sv, A_n, A_c, A_j, A_p, R_d)
      !   end do
      ! case ("de Boeck")
      !   ! Estimate Tleaf_C according to de Boeck (2012), iteratively
      !   ! solved with A_n and g_sto, starting with Tleaf_C = Tair_C.
      !   Tleaf_C = Tair_C
      !   call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                g_sv, A_n, A_c, A_j, A_p, R_d)
      !   do i = 1, MAX_ITERATIONS
      !     ! TODO: real "hypostomatous" setting, real "cloud cover" value?
      !     Tleaf_C = leaf_temp_de_Boeck(R, eact*1e3, Tair_C, Tleaf_C, P*1e3, &
      !                                  u, g_sv*1e-6, .true., Lm, albedo, 1.0, &
      !                                  pgc%Tleaf_balance_threshold, &
      !                                  pgc%Tleaf_adjustment_factor, &
      !                                  pgc%Tleaf_max_iterations)
      !     call farquhar_photosynthesis(Tleaf_C, CO2, eact*1e3, PPFD, g_bv, &
      !                                  pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, &
      !                                  g_sv, A_n, A_c, A_j, A_p, R_d)
      !   end do
      case ("Ewert")
          Tleaf_C = Tair_C

          call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARsun*4.57, g_bv, &
                                       pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
                                       fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
                                       g_sv, A_n_sun, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)

          call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARshade*4.57, g_bv, &
                                       pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
                                       fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
                                       g_sv, A_n_shade, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)
          call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARshade*4.57, g_bv, &
                                       pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
                                       fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
                                       g_sv, A_n, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)
          LAIsun = sunlit_LAI(LAI, sinB)
          LAIshade = LAI - LAIsun
          Canopy_A_n = LAIsun * A_n_sun + LAIshade * A_n_shade
      end select


      ! =================================== Skipping Select and just running EWERT ================================== !
      ! =================================== ----------------------------------------- ================================== !
      ! =================================== ----------------------------------------- ================================== !

      ! Tleaf_C = Tair_C

      ! call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARsun*4.57, g_bv, &
      !                             pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
      !                             fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
      !                             g_sv, A_n_sun, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)

      ! call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARshade*4.57, g_bv, &
      !                             pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
      !                             fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
      !                             g_sv, A_n_shade, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)
      ! call farquhar_photosynthesis_2(pgc, season, Tleaf_C, CO2, eact*1e3, PARshade*4.57, g_bv, &
      !                             pgc%g_sto_0, pgc%m, V_cmax_25, J_max_25, D_0*1e3, O3, td, dd, hr, &
      !                             fO3_h_1_hist, fO3_d_hist, daily_thermal_temperatures, &
      !                             g_sv, A_n, A_c, A_j, A_p, O3up, O3up_acc, fO3_h, fO3_d, R_d)
      ! LAIsun = sunlit_LAI(LAI, sinB)
      ! LAIshade = LAI - LAIsun
      ! Canopy_A_n = LAIsun * A_n_sun + LAIshade * A_n_shade

      ! =================================== ----------------------------------------- ================================== !

      ! Convert g_sto from umol to mmol, and from H2O to O3
      g_sto = DRATIO * (max(0.0, g_sv) / 1000)
    end subroutine gsto_pn

    ! !> Use quartic solution for leaf energy balance equation, from Nikolov (1995),
    ! !! to estimate the leaf temperature from ambient conditions.
    ! pure function nikolov_leaf_temperature(Tair_C, Pr, R_i, e_a, g_sv, g_bv, Wstat) result(Tleaf_C)
    !   real, intent(in) :: Tair_C    !< Ambient air temperature (degrees C)
    !   real, intent(in) :: Pr        !< Air pressure (Pa)
    !   real, intent(in) :: R_i       !< Net radiation absorbed by leaf (W m-2)
    !   real, intent(in) :: e_a       !< Water vapour pressure in ambient air (Pa)
    !   real, intent(in) :: g_sv      !< Leaf stomatal conductance to water vapour (micromol m-2 s-1)
    !   real, intent(in) :: g_bv      !< Leaf boundary layer conductance to water (micromol m-2 s-1)
    !   real, intent(in) :: Wstat     !< Leaf wetness status (0 = dry, 1 = wet)

    !   real :: Tleaf_C  !< Output: leaf temperature (degrees C)

    !   !> Specific heat capacity of dry air at standard pressure and 20C (J kg-1 K-1)
    !   real, parameter :: C_P = 1010.0
    !   !> Triple-point temperature of water (degrees K)
    !   real, parameter :: T3_H2O = 273.16
    !   !> Leaf thermal emissivity
    !   real, parameter :: LTE = 0.975
    !   !> Stefan-Boltzmann constant (W m-2 K-4)
    !   real, parameter :: SBC = 5.670373e-8

    !   ! Constants for saturation vapour pressure curve:
    !   !   e_s(T) = a*T**4 + b*T**3 + c*T**2 + d*T + e
    !   real, parameter :: SVP_A = 5.82436e-4
    !   real, parameter :: SVP_B = 1.5842e-2
    !   real, parameter :: SVP_C = 1.55186
    !   real, parameter :: SVP_D = 44.513596
    !   real, parameter :: SVP_E = 607.919

    !   ! For deriving the quartic coefficients
    !   real :: Cfm, lambda, psychro, Tvir, rho, h_e, h_t, k, a, b, c, d
    !   ! For solving the quartic
    !   real :: y, E, sa, P, Q, Dscr, R, t1

    !   ! Conversion from micromol m-2 s-1 to m s-1
    !   Cfm = 8.3089764 * 1e-6 * ((Tair_C + T0) / Pr)
    !   ! Latent heat of vapourisation for water (J kg-1)
    !   lambda = (-0.0000614342*Tair_C**3 + 0.00158927*Tair_C**2 - 2.36418*Tair_C + 2500.79) * 1000
    !   ! Psychrometric parameter (Pa C-1)
    !   psychro = (C_P * Pr) / (0.622 * lambda)
    !   ! Virtual temperature for density calcualation (K)
    !   Tvir = (Tair_C + T0) / (1 - (0.378 * (e_a / Pr)))
    !   ! Density of air (kg m-3)
    !   rho = Pr / (287.058 * Tvir)

    !   if (Wstat <= 0) then
    !     h_e = ((rho * C_P) / psychro) * Cfm * g_bv
    !   else
    !     h_e = ((rho * C_P) / psychro) * Cfm * ((g_sv * g_bv) / (g_sv + g_bv))
    !   end if

    !   ! 0.924 = ratio of heat conductance to vapour conductance
    !   h_t = rho * C_P * Cfm * (0.924 * g_bv)

    !   ! Set up coefficients for quartic: T**4 + a*T**3 + b*T**2 + c*T + d = 0
    !   k = 1 / ((2 * LTE * SBC) + (SVP_A * h_e))
    !   a = ((8  * LTE * SBC * T3_H2O   ) + (SVP_B * h_e)) * k
    !   b = ((12 * LTE * SBC * T3_H2O**2) + (SVP_C * h_e)) * k
    !   c = ((8  * LTE * SBC * T3_H2O**3) + (SVP_D * h_e) + h_t) * k
    !   d = ((2  * LTE * SBC * T3_H2O**4) + (SVP_E * h_e) - (h_t * Tair_C) - (h_e * e_a) - R_i) * k

    !   ! Presumably the analytical solution for the quartic?
    !   y = a * c - 4 * d
    !   E = b**2
    !   sa = a**2
    !   P = (3 * y - E) / 9
    !   Q = (b * (2 * E - 9 * y) - 27 * (d * (4 * b - sa) - c**2)) / 54
    !   Dscr = sqrt(Q**2 + P**3)
    !   y = exp(log(Q + Dscr) / 3) - exp(log(Dscr - Q) / 3) + b / 3
    !   R = sqrt(0.25 * sa + y - b)
    !   t1 = 0.25 * (a * (4 * b - sa) - 8 * c) / R
    !   E = 0.5 * sa - b - y
    !   Tleaf_C = -0.25 * a - 0.5 * (R - sqrt(E - t1))
    ! end function nikolov_leaf_temperature

  end module DO3SE_Photosynthesis_ml
