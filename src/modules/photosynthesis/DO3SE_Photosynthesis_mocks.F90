! saturated_vapour_pressure

module DO3SE_Photosynthesis_mocks

    implicit none

    public

    contains
      pure real function saturated_vapour_pressure(Ts_C)
        real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)

        saturated_vapour_pressure = 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))

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
