module mod_cplvar_sea
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public truncate_cplvar_sea
    public vsea_input, vsea_output

    ! Input and output sea variables exchanged by coupler
    ! Ocean model input variables
    type(rpe_var) :: vsea_input(ix*il,8)

    ! Ocean model output variablesend module
    type(rpe_var) :: vsea_output(ix*il,3)

    contains
        subroutine truncate_cplvar_sea()
            vsea_input = vsea_input
            vsea_output = vsea_output
        end subroutine
end module
