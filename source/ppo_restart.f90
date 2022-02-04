subroutine restart(jday)
    !  subroutine restart (jday)
    !
    !  Purpose : read or write a restart file
    !  Input :   JDAY  = 0 : read model variables from a restart file
    !                  > 0 : write model variables  to a restart file
    !                        at selected dates and at the end of run
    !

    use mod_atparam
    use mod_dynvar, only: vor, div, t, ps, tr
    use mod_date, only: iyear, imonth, iday, ihour
    use mod_downscaling, only: mx_in, nx_in, kx_in, ix_in, il_in, &
            calc_grid_weights
    use mod_prec
    use humidity,only: zero_c
    use rp_emulator

    implicit none

    !---------------------------------------------------------------------------
    integer, intent(in) :: jday
    integer :: yyyy, mm, dd, hh
    character(len=14) :: filename='yyyymmddhh.rst'

    ! Prognostic spectral variables. Can be at different resolution.
    ! Vorticity
    complex(dp) :: vor_in(mx_in, nx_in, kx_in, 2)

    ! Divergence
    complex(dp) :: div_in(mx_in, nx_in, kx_in, 2)

    ! Absolute temperature
    complex(dp) :: T_in(mx_in, nx_in, kx_in, 2)

    ! Log of (norm.) sfc pressure (p_s/p0)
    complex(dp) :: Ps_in(mx_in, nx_in, 2)

    ! Tracers (tr.1: spec. humidity in g/kg)
    complex(dp) :: tr_in(mx_in, nx_in, kx_in, 2, ntr)

    ! Truncation scale for initial conditions
    integer :: mx_tr, nx_tr

    !---------------------------------------------------------------------------
    if (jday==0) then
        ! 1. Read the restart dataset corresponding to the specified initial date
        write (filename(1:4),'(i4.4)') iyear
        write (filename(5:6),'(i2.2)') imonth
        write (filename(7:8),'(i2.2)') iday
        write (filename(9:10),'(i2.2)') ihour
        open (3, file=filename, form='unformatted')

        read (3,end=200) yyyy, mm, dd, hh

        print *, filename, yyyy,mm,dd,hh
        STOP

        if (yyyy/=iyear .or. mm/=imonth .or. dd/=iday .or. hh/=ihour) then
            print *, iyear, imonth, iday, ihour
            print *, yyyy, mm, dd, hh
            print*, "WARNING: Date in restart file does not match requested date"
            ! stop 'Date in restart file does not match requested date'
        end if

        print '(A,I4.4,A,I2.2,A,I2.2,A,I2.2)',&
                'Read restart dataset for year/month/date/hour: ', &
                iyear,'/',imonth,'/',iday,'/',ihour

        ! Load data at input resolution
        read (3) vor_in
        read (3) div_in
        read (3) T_in
        read (3) Ps_in
        read (3) tr_in

        ! If input vertical levels are different, interpolate
        if (kx/=kx_in) then
            ! TODO
            stop 'Input must have same vertical levels as output'
        end if

        ! Copy prognostic variables matching the truncation scale of input and
        ! model run
        mx_tr = min(mx, mx_in)
        nx_tr = min(nx, nx_in)

        vor = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        div = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        T   = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        Ps  = CMPLX(0.0_dp, 0.0_dp, kind=dp)
        tr  = CMPLX(0.0_dp, 0.0_dp, kind=dp)

        ! Reduce precision of input fields
        vor(1:mx_tr, 1:nx_tr, :, :)    = vor_in(1:mx_tr, 1:nx_tr, :, :)
        div(1:mx_tr, 1:nx_tr, :, :)    = div_in(1:mx_tr, 1:nx_tr, :, :)
        T  (1:mx_tr, 1:nx_tr, :, :)    = T_in  (1:mx_tr, 1:nx_tr, :, :)
        Ps (1:mx_tr, 1:nx_tr, :)       = Ps_in (1:mx_tr, 1:nx_tr, :)
        tr (1:mx_tr, 1:nx_tr, :, :, :) = tr_in (1:mx_tr, 1:nx_tr, :, :, :)

        ! Initialise gridded surface fields
        if (ix_in/=ix .or. il_in/=il) call calc_grid_weights()
        call rest_land(0)
        call rest_sea(0)
        close (3)
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
        vor_in(1:mx, 1:nx, :, :)    = vor(1:mx, 1:nx, :, :)
        vor_in = vor_in/ 3600.0_dp
        div_in(1:mx, 1:nx, :, :)    = div(1:mx, 1:nx, :, :)
        div_in = div_in / 3600.0_dp
        T_in  (1:mx, 1:nx, :, :)    = T  (1:mx, 1:nx, :, :)
        T_in(1,1,:,:) = T_in(1,1,:,:) + cmplx(sqrt(2.0_dp)*zero_c, kind=dp)

        Ps_in (1:mx, 1:nx, :)       = Ps (1:mx, 1:nx, :)
        tr_in (1:mx, 1:nx, :, :, :) = tr (1:mx, 1:nx, :, :, :)
        ! print*, tr_in(:,:,1,1,1)
        ! Write prognostic variables to restart file
        write (10) vor_in
        write (10) div_in
        write (10) t_in
        write (10) ps_in
        write (10) tr_in

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
end subroutine restart
