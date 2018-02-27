!> @brief
!> Constants for initialization of dynamics.
module mod_dyncon0
    implicit none

    namelist /dynamics/ gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    ! Ref. temperature lapse rate (-dT/dz in deg/km)
    real :: gamma = 6.0

    ! Ref. scale height for pressure (in km)
    real :: hscale = 7.5

    ! Ref. scale height for spec. humidity (in km)
    real :: hshum = 2.5

    ! Ref. relative humidity of near-surface air
    real :: refrh1 = 0.7

    ! Max damping time (in hours) for hor. diffusion (del^6) of temperature and
    ! vorticity
    real :: thd = 2.4

    ! Max damping time (in hours) for hor. diffusion (del^6)
    ! of divergence
    real :: thdd = 2.4

    ! Max damping time (in hours) for extra diffusion (del^2)
    ! in the stratosphere 
    real :: thds = 12.0

    ! Damping time (in hours) for drag on zonal-mean wind
    ! in the stratosphere 
    real :: tdrs = 24.0 * 30.0

    contains
        subroutine setup_dynamics(fid)
            integer, intent(in) :: fid

            read(fid, dynamics)

            write(*, dynamics)
        end subroutine setup_dynamics
end module
