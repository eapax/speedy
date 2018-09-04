module mod_fft
    use mod_atparam, only: ix
    use mod_prec, only: dp

    implicit none

    real(dp), allocatable :: wsave(:)

    contains
        subroutine setup_fft()
            allocate(wsave(2*ix+15))
        end subroutine setup_fft
end module
