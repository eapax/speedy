!> @brief
!> Constants for initialization of dynamics.
module mod_dyncon0
    use rp_emulator
    use mod_prec

    implicit none

    namelist /dynamics/ gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    ! Namelist parameters used to set up model constants
    ! Ref. temperature lapse rate (-dT/dz in deg/km)
    real(dp) :: gamma
    ! Ref. scale height for pressure (in km)
    real(dp) :: hscale
    ! Ref. scale height for spec. humidity (in km)
    real(dp) :: hshum
    ! Ref. relative humidity of near-surface air
    real(dp) :: refrh1
    ! Max damping time (in hours) for hor. diffusion (del^6) of temperature and
    ! vorticity
    real(dp) :: thd
    ! Max damping time (in hours) for hor. diffusion (del^6)
    ! of divergence
    real(dp) :: thdd
    ! Max damping time (in hours) for extra diffusion (del^2)
    ! in the stratosphere 
    real(dp) :: thds

    ! Namelist parameters used for model integration
    ! Damping time (in hours) for drag on zonal-mean wind
    ! in the stratosphere 
    type(rpe_var) :: tdrs

    contains
        subroutine setup_dynamics(fid)
            integer, intent(in) :: fid

            read(fid, dynamics)

            write(*, dynamics)
        end subroutine setup_dynamics
    
        subroutine truncate_dyncon0()
            call apply_truncation(tdrs)
        end subroutine truncate_dyncon0
end module
