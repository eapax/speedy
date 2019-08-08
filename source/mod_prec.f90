module mod_prec
    use rp_emulator
    use, intrinsic :: iso_fortran_env

    implicit none

    private
    public dp, sp, setup_precision, set_precision

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32
    integer :: &
            rp_convection=52, &
            rp_condensation=52, &
            rp_cloud=52, &
            rp_sw_radiation=52, &
            rp_lw_radiation=52, &
            rp_surface_fluxes=52, &
            rp_vertical_diffusion=52

    namelist /precisions/ &
            RPE_ACTIVE, RPE_IEEE_HALF, &
            rp_convection, rp_condensation, rp_cloud, &
            rp_sw_radiation, rp_lw_radiation, rp_surface_fluxes, &
            rp_vertical_diffusion

    ! Track previous precision
    integer :: rp_previous = 52

    contains

        subroutine setup_precision()
            ! Load values for precision in different parts of the model from
            ! a text file. This way I don't need to recompile the model every
            ! time I want to run with a different precision.
            open(99, file='precisions.nml')
            read(99, precisions)
            close(99)

            call set_precision('Double')
        end subroutine

        subroutine set_precision(mode)
            ! Set the global precision 'RPE_DEFAULT_SBITS' for specific parts of
            ! the model.
            character (len=*), intent(in) :: mode

            ! Save current precision before switching
            rp_previous = RPE_DEFAULT_SBITS

            select case(mode)
                case('Double')
                RPE_DEFAULT_SBITS = 52

                case('Single')
                RPE_DEFAULT_SBITS = 23

                case('Half')
                RPE_DEFAULT_SBITS = 10

                case('Convection')
                RPE_DEFAULT_SBITS = rp_convection

                case('Condensation')
                RPE_DEFAULT_SBITS = rp_condensation

                case('Cloud')
                RPE_DEFAULT_SBITS = rp_cloud

                case('Short-Wave Radiation')
                RPE_DEFAULT_SBITS = rp_sw_radiation

                case('Long-Wave Radiation')
                RPE_DEFAULT_SBITS = rp_lw_radiation

                case('Surface Fluxes')
                RPE_DEFAULT_SBITS = rp_surface_fluxes

                case('Vertical Diffusion')
                RPE_DEFAULT_SBITS = rp_vertical_diffusion

                case('Previous')
                RPE_DEFAULT_SBITS = rp_previous

                case Default
                    print *,'Trying to use an unset precision ', mode
                    stop
            end select
        end subroutine
end module
