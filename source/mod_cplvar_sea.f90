module mod_cplvar_sea
    use mod_atparam

    implicit none

    ! Input and output sea variables exchanged by coupler
    ! Ocean model input variables
    real, allocatable :: vsea_input(:,:)

    ! Ocean model output variablesend module
    real, allocatable :: vsea_output(:,:)
    
    contains
    
        subroutine setup_cplvar_sea()
            allocate(vsea_input(ix*il,8))
            allocate(vsea_output(ix*il,3))
        end subroutine setup_cplvar_sea
end module
