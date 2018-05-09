module mod_cplvar_sea
    use mod_atparam
    use rp_emulator

    implicit none

    ! Input and output sea variables exchanged by coupler
    ! Ocean model input variables
    type(rpe_var), allocatable :: vsea_input(:,:)

    ! Ocean model output variablesend module
    type(rpe_var), allocatable :: vsea_output(:,:)

    contains

        subroutine setup_cplvar_sea()
            allocate(vsea_input(ngp,8))
            allocate(vsea_output(ngp,3))
        end subroutine setup_cplvar_sea

        subroutine truncate_cplvar_sea()
            vsea_input = vsea_input
            vsea_output = vsea_output
        end subroutine truncate_cplvar_sea
end module
