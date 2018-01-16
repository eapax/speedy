module mod_prec
    use rp_emulator
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public dp, sp, setup_precision, set_precision

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32
    integer :: reduced_precision, rp_initial_values, rp_spectral_transform, &
            rp_grid_physics, rp_grid_dynamics, rp_spectral_dynamics, &
            rp_diffusion, rp_timestepping, rp_prognostics, rp_tendencies, &
            rp_initialisation, rp_parameters

    contains

        subroutine setup_precision()
            ! Load values for precision in different parts of the model from
            ! a text file. This way I don't need to recompile the model every
            ! time I want to run with a different precision.

            open(99, file='precision.txt')
            read (99,*) reduced_precision
            read (99,*) rp_initial_values
            read (99,*) rp_spectral_transform
            read (99,*) rp_grid_physics
            read (99,*) rp_grid_dynamics
            read (99,*) rp_spectral_dynamics
            read (99,*) rp_diffusion
            read (99,*) rp_timestepping
            read (99,*) rp_prognostics
            read (99,*) rp_tendencies
            read (99,*) rp_initialisation
            read (99,*) rp_parameters
            close(99)

            call set_precision('Default')
        end subroutine

        subroutine set_precision(mode)
            ! Set the global precision 'RPE_DEFAULT_SBITS' for specific parts of
            ! the model.
            character (len=*), intent(in) :: mode

            select case(mode)
                case('Default')
                RPE_DEFAULT_SBITS = reduced_precision

                case('Full')
                RPE_DEFAULT_SBITS = 52

                case('Single')
                RPE_DEFAULT_SBITS = 23

                case('Half')
                RPE_DEFAULT_SBITS = 10

                case('Initial Values')
                RPE_DEFAULT_SBITS = rp_initial_values

                case('Spectral Transform')
                RPE_DEFAULT_SBITS = rp_spectral_transform

                case('Grid Physics')
                RPE_DEFAULT_SBITS = rp_grid_physics

                case('Grid Dynamics')
                RPE_DEFAULT_SBITS = rp_grid_dynamics

                case('Spectral Dynamics')
                RPE_DEFAULT_SBITS = rp_spectral_dynamics

                case('Diffusion')
                RPE_DEFAULT_SBITS = rp_diffusion

                case('Timestepping')
                RPE_DEFAULT_SBITS = rp_timestepping

                case('Prognostics')
                RPE_DEFAULT_SBITS = rp_prognostics

                case('Tendencies')
                RPE_DEFAULT_SBITS = rp_tendencies

                case('Initialisation')
                RPE_DEFAULT_SBITS = rp_initialisation

                case('Parameters')
                RPE_DEFAULT_SBITS = rp_parameters
            end select
        end subroutine
end module
