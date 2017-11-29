module mod_prec
    use mod_atparam, only: ix, il, kx
    use rp_emulator
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public dp, sp, setup_precision, set_precision, set_precision_grid, &
            set_precision_spectral, determine_location_dependent_precision

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32
    integer :: reduced_precision
    integer :: zeroth_mode_precision
    integer :: grid_dynamics_precision
    integer :: location_dependent_precision(ix, il, kx)
    integer :: initial_precision

    contains

        subroutine setup_precision()
            ! Load values for precision in different parts of the model from
            ! a text file. This way I don't need to recompile the model every
            ! time I want to run with a different precision.

            open(99, file='precision.txt')
            read (99,*) reduced_precision
            read (99,*) zeroth_mode_precision
            read (99,*) grid_dynamics_precision
            read (99,*) initial_precision

            print *, 'Reduced precision = ', reduced_precision
            print *, 'Zeroth mode precision = ', zeroth_mode_precision
            print *, 'Grid-dynamics precision = ', grid_dynamics_precision
            print *, 'Initial precision = ', initial_precision

            location_dependent_precision = grid_dynamics_precision

            call set_precision(0)
        end subroutine

        subroutine determine_location_dependent_precision(tendency, threshold)
            ! Determing the precision to use at each gridpoint based on the
            ! magnitude of the physics tendencies
            type(rpe_var), intent(in) :: tendency(ix, il, kx)
            real, intent(in) :: threshold
            integer :: i, j, k, n

            n=0

            do k=1,kx
                do j=1,il
                    do i=1,ix
                        if (abs(tendency(i,j,k)) > threshold) then
                            location_dependent_precision(i,j,k) = &
                                    grid_dynamics_precision
                            n=n+1
                        else
                            location_dependent_precision(i,j,k) = &
                                    reduced_precision
                        end if
                    end do
                end do
            end do

            print *, 'Number of reduced precision gridpoints = ', n

        end subroutine

        subroutine set_precision(n)
            ! Set the global precision 'RPE_DEFAULT_SBITS' for specific parts of
            ! the model.
            integer, intent(in) :: n

            select case(n)
                case default
                RPE_DEFAULT_SBITS = reduced_precision

                case(1)
                RPE_DEFAULT_SBITS = grid_dynamics_precision

                case(2)
                RPE_DEFAULT_SBITS = initial_precision
            end select
        end subroutine

        subroutine set_precision_grid(i, j, k)
            ! Set the global precision 'RPE_DEFAULT_SBITS' within a loop over
            ! gridpoints (i, j).
            integer, intent(in) :: i, j, k

            RPE_DEFAULT_SBITS = location_dependent_precision(i,j,k)
        end subroutine

        subroutine set_precision_spectral(m, n)
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
