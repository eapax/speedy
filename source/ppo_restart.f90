subroutine restart(jday)
    !  subroutine restart (jday)
    !
    !  Purpose : read or write a restart file
    !  Input :   JDAY  = 0 : read model variables from a restart file
    !                  > 0 : write model variables  to a restart file
    !                        at selected dates and at the end of run 
    !
    
    use mod_tsteps, only: nmonrs, iyear0, imont0
    use mod_atparam
    use mod_dynvar
    use mod_date, only: iyear, imonth, iday, ndaytot, ihour
    use rp_emulator
    use mod_prec, only: set_precision

    implicit none

    integer, intent(in) :: jday
    integer :: yyyy, mm, dd, hh, m, n
    character(len=14) :: filename='yyyymmddhh.rst'

    if (jday.eq.0) then
        ! 1. Read the restart dataset corresponding to the specified initial date
        write (filename(1:4),'(i4.4)') iyear
        write (filename(5:6),'(i2.2)') imonth
        write (filename(7:8),'(i2.2)') iday
        write (filename(9:10),'(i2.2)') ihour
        open (3, file=filename, form='unformatted')

        read (3,end=200) yyyy, mm, dd, hh

        if (yyyy/=iyear .or. mm/=imonth .or. dd/=iday .or. hh/=ihour) then
            stop 'Date in restart file does not match requested date'
        end if

        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
                'Read restart dataset for year/month/date/hour: ', &
                iyear,'/',imonth,'/',iday,'/',ihour

        ! Load data in full precision
        call set_precision('Full')

        read (3) vor
        read (3) div
        read (3) t
        read (3) ps
        read (3) tr

        ! Reduce precision of input fields
        call set_precision('Initial Values')
        vor = vor
        div = div
        t = t
        ps = ps
        tr = tr

        call rest_land(0)
        call rest_sea(0)
        close (3)

        call set_precision('Initialisation')
    else
        ! 2. Write date and model variables to the restart file
        print*, 'Write restart dataset for year/month/date/hour: ', &
                iyear,'/',imonth,'/',iday,'/',ihour

        ! Set filename to restart date
        write (filename(1:4),'(i4.4)') iyear
        write (filename(5:6),'(i2.2)') imonth
        write (filename(7:8),'(i2.2)') iday
        write (filename(9:10),'(i2.2)') ihour
        open (10, file=filename, form='unformatted')

        ! Write date to restart file
        write (10) iyear, imonth, iday, ihour

        ! Write prognostic variables to restart file
        write (10) vor
        write (10) div
        write (10) t
        write (10) ps
        write (10) tr

        ! Write surface fields to restart file
        call rest_land(1)
        call rest_sea(1)

        close(10)
    end if

    return

    ! 4. Stop integration if restart file is not found
    200 continue

    print*, ' No restart dataset for the specified initial date'

    stop 'invalid restart'
end
