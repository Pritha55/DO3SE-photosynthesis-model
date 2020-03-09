! saturated_vapour_pressure

module mocked_functions

  implicit none

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
end module mocked_functions
