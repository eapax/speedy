program agcm_main
    use mod_tsteps, only: nmonrs
    use mod_date, only: imonth, iday, ndays
    use rp_emulator
    use mod_prec, only: setup_precision

    implicit none

    ! program : agcm_main

    integer :: jday

    ! Setup reduced precision emulator
    RPE_IEEE_HALF = .false.
    RPE_ACTIVE = .true.
    call setup_precision()

    ! 1. initialization
    ! ndays = no. of integration days, set by agcm_init
    call agcm_init()

    print *, 'integration length in days: ', ndays

    ! 2. do loop over total no. of integration days
    do jday = 1, ndays
        ! 2.1 run atmospheric model for 1 day
        call agcm_1day(jday)

        ! 2.2 exchange data with coupler
        call agcm_to_coupler(jday)
        call coupler_to_agcm(jday)
    enddo

    ! Write restart file at end of run if not already written
    if (mod(imonth, nmonrs) /= 0 .or. iday /= 1) call restart(2)
end

subroutine agcm_1day(jday)
    ! subroutine agcm_1day (jday)
    !
    ! perform atm. model integration for 1 day, 
    ! post-proc. and i/o at selected times 

    use mod_tsteps, only: nsteps, nmonrs
    use mod_date, only: iyear, imonth, iday

    implicit none

    integer, intent(in) :: jday
    integer :: istep

    if (iday == 1) print *, ' start of year/month = ', iyear, imonth

    istep = 1 + (jday - 1) * nsteps

    ! 1. set forcing terms according to date
    call fordate(1)

    ! 2. set daily-average flux arrays to zero
    call dmflux(0)

    ! 3. integrate the atmospheric model for 1 day
    call stloop(istep)

    ! 4. Write restart file at the end of selected months
    if (iday == 1) then
        if (mod(imonth, nmonrs) == 0) call restart(2)
    endif
end
