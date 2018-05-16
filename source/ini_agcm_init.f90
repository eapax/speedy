subroutine agcm_init(ndays)
    !   purpose: initialization of atmos. model and coupling interface 
    !

    use mod_cpl_flags, only: icsea, isstan
    use mod_tsteps
    use mod_date, only: newdate, ndaytot, iyear, imonth, iday, ihour
    use ppo_output_stream, only: initialise_output
    use mod_prec, only: set_precision

    implicit none

    ! output:
    integer, intent(inout) :: ndays       ! total no. of integration days

    print *, ' hallo from speedy_agcm'

    ! 0. Read namelist and allocate arrays
    call ini_namelist()

    ! Initialise reduced precision constants
    call set_precision('Initialisation')
    call ini_rp()

    ! 1. set run initial time, duration, time-stepping and coupling options
    read (2,*) istart

    ! Read date from fort.2 file
    read (2,*) iyear
    read (2,*) imonth
    read (2,*) iday
    read (2,*) ihour
    iyear0 = iyear
    imont0 = imonth

    call newdate(0)

    print *, 'start date ', iyear, imonth, iday, ihour

    isst0 = (iyear - issty0) * 12 + imonth

    ndays = ndaytot

    ! check consistency of coupling and prescribed SST anomaly flags
    if (icsea >= 4) isstan = 1

    ! 2. initialization of atmospheric model constants and variables
    call ini_atm()

    ! 3. initialization of coupled modules (land, sea, ice)
    call ini_coupler(istart)

    ! 4. set up the forcing fields for the first time step
    call fordate(0)

    ! 5. do the initial (2nd-order) time step, initialize the semi-impl. scheme
    call stepone()

    ! 6. Set up model output
    call initialise_output()

    ! Truncate parameters and derived constants
    call set_precision('Parameters')
    call truncate_rp()
    call set_precision('Default')
end subroutine
