!> @brief
!> Constants for initialization of dynamics.
module mod_dyncon0
    use mod_prec

    implicit none

    private
    public gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    ! Ref. temperature lapse rate (-dT/dz in deg/km)
    real(dp), parameter :: gamma = 6.0

    ! Ref. scale height for pressure (in km)
    real(dp), parameter :: hscale = 7.5

    ! Ref. scale height for spec. humidity (in km)
    real(dp), parameter :: hshum = 2.5

    ! Ref. relative humidity of near-surface air
    real(dp), parameter :: refrh1 = 0.7

    ! Max damping time (in hours) for hor. diffusion (del^6) of temperature and
    ! vorticity
    real(dp), parameter :: thd = 2.4

    ! Max damping time (in hours) for hor. diffusion (del^6)
    ! of divergence
    real(dp), parameter :: thdd = 2.4

    ! Max damping time (in hours) for extra diffusion (del^2)
    ! in the stratosphere 
    real(dp), parameter :: thds = 12.0

    ! Damping time (in hours) for drag on zonal-mean wind
    ! in the stratosphere 
    real(dp), parameter :: tdrs = 24.0 * 30.0
end module
