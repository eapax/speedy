subroutine agcm_init()
    !   purpose: initialization of atmos. model and coupling interface

    use mod_cpl_flags, only: icsea, isstan
    use mod_date, only: ini_date, istart
    use ppo_output_stream, only: initialise_output
    use mod_fordate, only: ini_fordate

    implicit none

    print *, ' hallo from speedy_agcm'

    ! 0. Read namelist and allocate arrays
    call ini_namelist()

    ! 1. set run initial time, duration, time-stepping and coupling options
    call ini_date()

    ! check consistency of coupling and prescribed SST anomaly flags
    if (icsea>=4) isstan = 1

    ! 2. initialization of atmospheric model constants and variables
    call ini_atm()

    ! 3. initialization of coupled modules (land, sea, ice)
    call ini_coupler(istart)

    ! 4. set up the forcing fields for the first time step
    call ini_fordate()

    ! 5. do the initial (2nd-order) time step, initialize the semi-impl. scheme
    call stepone()

    ! 6. Set up model output
    call initialise_output()
end subroutine
