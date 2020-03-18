module DO3SE_Photosynthesis_helpers_ml

    implicit none
    private

    public :: calc_fO3_h
    public :: calc_f_LA
    public :: calc_f_LS
    public :: temp_dep
    public :: temp_dep_inhibit

contains
  pure function calc_fO3_h(O3up) result(fO3_h)
    real, intent(in) :: O3up
    real :: fO3_h
    real, parameter :: Gamma_1 = 0.060
    real, parameter :: Gamma_2 = 0.0045

    if (O3up <= (Gamma_1/Gamma_2)) then
      fO3_h = 1
    else if ((Gamma_1/Gamma_2) < O3up .AND. O3up < ((1 + Gamma_1) / Gamma_2)) then
      fO3_h = 1 + Gamma_1 - Gamma_2*O3up
    else if (O3up >= ((1 + Gamma_1) / Gamma_2)) then
      fO3_h = 0
    end if
  end function calc_fO3_h

  pure function calc_f_LA(t_lem, t_lma, t_l, td) result(f_LA)
    real, intent(in) :: t_lem
    real, intent(in) :: t_lma
    real, intent(in) :: t_l
    real, intent(in) :: td
    real :: f_LA

    if (td < t_lem) then
      f_LA = 1
    else if (t_lem <= td .AND. td < t_l) then
      f_LA = 1 - (td - t_lem) / t_lma
    else if (td >= t_l) then
      f_LA = 0
    end if
  end function

  pure function calc_f_LS(t_lem, t_lep, t_lma, t_l, fO3_l, td) result(f_LS)
    real, intent(in) :: t_lem
    real, intent(in) :: t_lep
    real, intent(in) :: t_lma
    real, intent(in) :: t_l
    real, intent(in) :: fO3_l
    real, intent(in) :: td
    real :: f_LS

    if (td <= (t_lem + t_lep)) then
      f_LS = 1
    else if ((t_lem + t_lep) < td .AND. td < t_l) then
      f_LS = 1 - (td - t_lem - t_lep) / (t_lma/fO3_l - t_lep)
    else if (td > t_l) then
      f_LS = 0
    end if
  end function


    !> Calculate parameter from temperature dependence curve
  pure function temp_dep(P_ref, T_ref, H_a, T, R) result (P)
    real, intent(in) :: P_ref   ! Parameter value at T_ref
    real, intent(in) :: T_ref   ! Reference temperature (degrees K)
    real, intent(in) :: H_a     ! Activation energy (J/mol)
    real, intent(in) :: T       ! Temperature (degrees K)
    real, intent(in) :: R           !< Global radiation (W m-2)

    real :: P

    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T))
  end function temp_dep

  !> Calculate parameter from temperature dependence curve with high
  !! temperature inhibition
  pure function temp_dep_inhibit(P_ref, T_ref, H_a, H_d, S, T, R) result (P)
    real, intent(in) :: P_ref   ! Parameter value at T_ref
    real, intent(in) :: T_ref   ! Reference temperature (degrees K)
    real, intent(in) :: H_a     ! Activation energy (J/mol)
    real, intent(in) :: H_d     ! Deactivation energy (J/mol)
    real, intent(in) :: S       ! Entropy term (J/(mol*K))
    real, intent(in) :: T       ! Temperature (degrees K)
    real, intent(in) :: R           !< Global radiation (W m-2)


    real :: P

    P = P_ref * exp((H_a * (T - T_ref)) / (T_ref * R * T)) &
        * (  (1 + exp((T_ref*S - H_d) / (T_ref*R))) &
          / (1 + exp((T*S - H_d) / (T*R))) )
  end function temp_dep_inhibit

end module DO3SE_Photosynthesis_helpers_ml
