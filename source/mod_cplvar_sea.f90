module mod_cplvar_sea
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Input and output sea variables exchanged by coupler
    ! Ocean model input variables
    real(dp), allocatable :: vsea_input(:,:)

    ! Ocean model output variablesend module
    real(dp), allocatable :: vsea_output(:,:)

    contains

        subroutine setup_cplvar_sea()
            allocate(vsea_input(ngp,8))
            allocate(vsea_output(ngp,3))
        end subroutine setup_cplvar_sea
end module
