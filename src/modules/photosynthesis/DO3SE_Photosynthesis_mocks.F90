! saturated_vapour_pressure

module DO3SE_Photosynthesis_mocks

    implicit none

    public
    ! ! PHYSICAL CONSTANTS
    ! ! real, parameter :: PI = 3.141592653589793238
    ! ! real, parameter :: DEG2RAD = 0.017453292519943295

    ! ! !> Atmospheric pressure at sea level (kPa)
    ! ! real, parameter :: seaP = 101.325

    ! !> 0 degrees Celsius in Kelvin
    ! real, parameter :: T0 = 273.15
    ! !> Stefan-Boltzmann constant (W m-2 K-4)
    ! ! real, parameter :: SBC = 5.670373e-8

    ! !> Universal gas constant (J K-1 mol-1)
    ! ! TODO: update this value to 8.3144621
    ! real, parameter :: R = 8.314472

    ! !> Approximate fraction of global radiation in PAR waveband
    ! ! real, parameter :: PARfrac = 0.45
    ! !> PAR conversion from W m-2 to umol photons m-2 s-1
    ! ! real, parameter :: PAR_Wm2_to_photons = 4.57
    ! !> Net radiation conversion from MJ m-2 h-1 to W m-2
    ! real, parameter :: Rn_MJ_to_W = 277.8

    ! !> Molecular diffusivity of O3 in air (m2 s-1)
    ! ! real, parameter :: DIFF_O3 = 0.000015
    ! !> Molecular diffusivity of H2O (m2 s-1)
    ! ! real, parameter :: DIFF_H2O = 0.000025
    ! !> Ratio between molecular diffusivity of O3 and H2O
    ! real, parameter :: DRATIO = 0.663


    ! ! MODEL CONSTANTS
    ! REAL, parameter :: UNDEF = -999.0
    ! ! INTEGER, parameter :: IUNDEF = -999

    ! ! TODO: see how much we can remove these
    ! ! INTEGER, parameter :: MAX_LC = 3      !< Maximum number of land covers (used in some static allocations)
    ! ! INTEGER, parameter :: MAX_LAYERS = 5  !< Maximum number of layers (used in some static allocations)

    ! ! real, parameter :: DT = 60*60   !< Number of seconds in a timestep

    ! ! !> Canopy displacement (fraction of canopy height)
    ! ! real, parameter :: CANOPY_D = 0.7
    ! ! !> Canopy roughness length (fraction of canopy height)
    ! ! real, parameter :: CANOPY_Z0 = 0.1

    ! ! !> ASW (available soil water) for minimum gsto (percent of ASW at field capacity)
    ! ! real, parameter :: ASW_MIN = 0.0
    ! ! !> ASW (available soil water) for maximum gsto (percent of ASW at field capacity)
    ! ! real, parameter :: ASW_MAX = 50.0

    ! ! real, parameter :: LEAF_G_HEAT = 0.135
    ! real, parameter :: LEAF_G_H2O = 0.147
    ! ! real, parameter :: LEAF_G_CO2 = 0.110
    ! ! real, parameter :: LEAF_G_O3 = 0.105


    contains
      pure function saturated_vapour_pressure(Tleaf_C) RESULT(j)
          real, intent(in) :: Tleaf_C ! input
        real             :: j ! output
        j = 1

      end function saturated_vapour_pressure

      pure function sunlit_LAI(LAI, sinB) RESULT(j)
        real, intent(in) :: LAI, sinB ! input
        real             :: j ! output
        j = 1

      end function sunlit_LAI

      ! pure subroutine gsto_pn(g_sto)
      !   real, intent(out) :: g_sto      !< Output: Stomatal conductance (mmol m-2 PLA s-1)
      !   g_sto = 1
      ! end subroutine gsto_pn
  end module DO3SE_Photosynthesis_mocks
