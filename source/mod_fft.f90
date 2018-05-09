module mod_fft
    use mod_atparam, only: ix
    use rp_emulator

    implicit none

    type(rpe_var), allocatable :: wsave(:)

    contains
        subroutine setup_fft()
            allocate(wsave(2*ix+15))
        end subroutine setup_fft

        subroutine truncate_fft()
            wsave(1:2*ix) = wsave(1:2*ix)
        end subroutine
end module
