module mod_prec
    use rp_emulator
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public dp, sp, reduced_precision, set_precision

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32
    integer, parameter :: zeroth_mode_precision = 16
    integer, parameter :: reduced_precision = 16

    contains

    subroutine set_precision(m, n)
        integer, intent(in) :: m, n

        if (m == 1 .and. n == 1) then
            RPE_DEFAULT_SBITS = zeroth_mode_precision
        else
            RPE_DEFAULT_SBITS = reduced_precision
        end if

    end subroutine
end module
