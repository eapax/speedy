module mod_fft
    use mod_atparam, only: ix
    use rp_emulator

    implicit none

    private
    public truncate_fft
    public wsave

    type(rpe_var) :: wsave(2*ix+15)

    contains
        subroutine truncate_fft()
            wsave = wsave
        end subroutine
end module
