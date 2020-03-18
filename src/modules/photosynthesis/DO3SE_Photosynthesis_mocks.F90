! saturated_vapour_pressure

module DO3SE_Photosynthesis_mocks

    implicit none

    public

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
