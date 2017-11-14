module mod_prec
    use rp_emulator
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public dp, sp, set_precision, setup_precision

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32
    integer :: zeroth_mode_precision
    integer :: reduced_precision

    contains

        subroutine setup_precision()
            ! Load values for precision in different parts of the model from
            ! a text file. This way I don't need to recompile the model every
            ! time I want to run with a different precision.

            open(99, file='precision.txt')
            read (99,*) reduced_precision
            read (99,*) zeroth_mode_precision

            call set_precision(0, 0)
        end subroutine

        subroutine set_precision(m, n)
            ! Set the global precision 'RPE_DEFAULT_SBITS' within a loop over
            ! wavenumbers (m, n). For the zeroth mode (m=1, n=1) the precision
            ! can be increased to a pre-defined precision. Otherwise use the
            ! reduced precision.
            integer, intent(in) :: m, n

            if (m == 1 .and. n == 1) then
                RPE_DEFAULT_SBITS = zeroth_mode_precision
            else
                RPE_DEFAULT_SBITS = reduced_precision
            end if

        end subroutine
end module
