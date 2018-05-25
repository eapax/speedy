!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    namelist /timestepping/ nsteps, nstdia, nmonrs, rob, wil

    ! No. of time steps in one day
    integer :: nsteps = 36

    ! Period (no. of steps) for diagnostic print-out
    integer :: nstdia = 9

    ! Period (no. of months) for restart file update
    integer :: nmonrs = 3

    ! Time step in seconds
    integer :: idelt
    type(rpe_var) :: delt
    
    ! 2 * time step in seconds
    type(rpe_var) :: delt2

    ! Damping factor in Robert time filter
    type(rpe_var) :: rob

    ! Parameter of Williams filter
    type(rpe_var) :: wil

    ! Coefficient for semi-implicit computations
    type(rpe_var) :: alph

    contains
        subroutine setup_timestepping(fid)
            integer, intent(in) :: fid

            read(fid, timestepping)
            idelt = 86400 / nsteps
            delt = real(idelt, dp)
            delt2 = 2.0 * idelt

            write(*, timestepping)
        end subroutine setup_timestepping

        subroutine truncate_tsteps()
            call apply_truncation(delt)
            call apply_truncation(delt2)
            call apply_truncation(rob)
            call apply_truncation(wil)
            call apply_truncation(alph)
        end subroutine truncate_tsteps
end module
