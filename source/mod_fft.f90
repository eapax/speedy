module mod_fft
    use mod_atparam, only: ix

    implicit none

    real, allocatable :: wsave(:)

    contains
        subroutine setup_fft()
            allocate(wsave(2*ix+15))
        end subroutine setup_fft
end module
