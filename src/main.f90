program DO3SE_main

    ! use DO3SE_Photosynthesis_ml

    implicit none

    character(len=256) :: arg1

    print *, "hello world"
    call get_command_argument(1, arg1)
    print * , arg1

    ! gsto_pn()
  
end program DO3SE_main
  