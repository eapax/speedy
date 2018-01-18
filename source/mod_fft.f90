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
            wsave(1:192) = wsave(1:192)
            wsave(196:) = wsave(196:)
        end subroutine
end module
