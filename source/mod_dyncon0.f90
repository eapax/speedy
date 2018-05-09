!> @brief
!> Constants for initialization of dynamics.
module mod_dyncon0
    use rp_emulator
    use mod_prec

    implicit none

    namelist /dynamics/ gamma, hscale, hshum, refrh1, thd, thdd, thds, tdrs

    ! Ref. temperature lapse rate (-dT/dz in deg/km)
    type(rpe_var) :: gamma

    ! Ref. scale height for pressure (in km)
    type(rpe_var) :: hscale

    ! Ref. scale height for spec. humidity (in km)
    type(rpe_var) :: hshum

    ! Ref. relative humidity of near-surface air
    type(rpe_var) :: refrh1

    ! Max damping time (in hours) for hor. diffusion (del^6) of temperature and
    ! vorticity
    type(rpe_var) :: thd

    ! Max damping time (in hours) for hor. diffusion (del^6)
    ! of divergence
    type(rpe_var) :: thdd

    ! Max damping time (in hours) for extra diffusion (del^2)
    ! in the stratosphere 
    type(rpe_var) :: thds

    ! Damping time (in hours) for drag on zonal-mean wind
    ! in the stratosphere 
    type(rpe_var) :: tdrs

    contains
        subroutine setup_dynamics(fid)
            integer, intent(in) :: fid

            read(fid, dynamics)

            write(*, dynamics)
        end subroutine setup_dynamics

        subroutine init_dyncon0
            gamma = gamma
            hscale = hscale
            hshum = hshum
            refrh1 = refrh1
            thd = thd
            thdd = thdd
            thds = thds
            tdrs = tdrs
        end subroutine
end module
