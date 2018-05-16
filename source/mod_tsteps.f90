!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
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
    real :: delt
    
    ! 2 * time step in seconds
    real :: delt2

    ! Damping factor in Robert time filter
    real :: rob = 0.05

    ! Parameter of Williams filter
    real :: wil = 0.53

    ! Coefficient for semi-implicit computations
    real :: alph

    contains
        subroutine setup_timestepping(fid)
            integer, intent(in) :: fid

            read(fid, timestepping)
            idelt = 86400 / nsteps
            delt = real(idelt)
            delt2 = 2.0 * idelt

            write(*, timestepping)
        end subroutine setup_timestepping
end module
