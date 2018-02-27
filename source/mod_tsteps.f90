!> @brief
!> Length of the integration and time stepping constants.
module mod_tsteps
    implicit none

    namelist /timestepping/ nmonts, ndaysl, nsteps, nstdia, nmonrs, iseasc,&
            rob, wil

    ! Integration length in months
    integer :: nmonts = 3

    ! No. of days in the last month of int. (max=30)
    integer :: ndaysl = 0

    ! No. of time steps in one day
    integer :: nsteps = 36

    ! Period (no. of steps) for diagnostic print-out
    integer :: nstdia = 9

    ! Period (no. of months) for restart file update
    integer :: nmonrs = 3

    ! Seasonal cycle flag (0=no, 1=yes)
    integer :: iseasc = 1

    ! Start flag (0: from rest, 1: from restart file)
    integer :: istart

    ! Year of initial date (4-digit, eg 1900)
    integer :: iyear0
    
    ! Month of initial date (1 to 12)
    integer :: imont0

    integer, parameter :: issty0 = 1979

    ! Record in SST anomaly file corr. to the initial month
    ! Initialized in agcm_init
    integer :: isst0
    
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
