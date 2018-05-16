module mod_date

    implicit none

    private
    public setup_date, ini_date, update_time, month_start
    public istart, iyear0, imonth0
    public iyear, imonth, iday, ihour, ndays
    public imont1, tmonth, tyear, issty0, isst0

    namelist /date/ istart, iyear0, imonth0, iday0, nmonts, ndaysl, iseasc

    ! Start flag (0: from rest, 1: from restart file)
    integer :: istart

    ! Initial date (read in Namelist)
    integer :: iyear0, imonth0, iday0

    ! Run length (read in Namelist)
    ! Integration length in months
    integer :: nmonts
    ! No. of days in the last month of int. (max=30)
    integer :: ndaysl

    ! Seasonal cycle flag (0=no, 1=yes)
    integer :: iseasc

    ! Run length in days (calculated in ini_date)
    integer :: ndays

    ! Current date variables (updated in newdate each day)
    integer :: iyear, imonth, iday
    ! Current time (update in update time each timestep)
    integer :: ihour=0, isecond=0

    ! Additional variables calculated in newdate to use for boundary forcing
    integer :: imont1
    real :: tmonth, tyear

    ! Record in SST anomaly file corr. to the initial month
    ! Initialized in ini_date
    integer, parameter :: issty0 = 1979
    integer :: isst0

    ! Calendar set-up (initialized in ini_date)
    integer :: ndaycal(12,2)

    ! 365-day calendar
    integer, parameter :: ncal = 365
    integer, parameter :: ncal365(12) = &
            (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

    contains
        subroutine setup_date(fid)
            integer, intent(in) :: fid

            read(fid, date)

            iyear=iyear0
            imonth=imonth0
            iday=iday0
        end subroutine setup_date

        subroutine ini_date()
            ! initilialize date variables
            integer :: jm, im

            ! calendar
            if (ncal == 365) then
                ndaycal(:,1) = ncal365(:)
            else
                ndaycal(:,1) = 30
            end if

            ! ndaycal(:,2) stores cumulative number of days up to that month
            ndaycal(1,2) = 0
            do jm = 2, 12
                ndaycal(jm,2) = ndaycal(jm-1,1)+ndaycal(jm-1,2)
            end do

            ! total no. of integration days
            ndays = ndaysl
            im = imonth0

            do jm=1,nmonts
                ndays = ndays+ndaycal(im,1)
                im = im+1
                if (im.eq.13) im=1
            end do

            print *, 'start date ', iyear, imonth, iday, ihour

            ! Find index of sst anomaly in file
            isst0 = (iyear - issty0) * 12 + imonth
        end subroutine

        subroutine update_time()
            ! Increment the time by a single timestep
            use mod_tsteps, only: idelt

            isecond = isecond + idelt
            do while(isecond >= 3600)
                ! Transfer one hour from the seconds counter to the hour counter
                isecond = isecond - 3600

                ! Check for a new day
                ihour = mod(ihour + 1, 24)
                if (ihour==0) call newdate()
            end do
        end subroutine update_time

        function month_start()
            ! Returns true if it is the first timestep in a month and false
            ! otherwise
            use mod_tsteps, only: idelt

            logical :: month_start

            month_start = (isecond < idelt .and. ihour == 0 .and. iday == 1)
        end function month_start

        subroutine newdate()
            ! set new date
            iday = iday+1

            ! Leap year and February?
            if (mod(iyear,4) == 0 .and. imonth == 2) then
                if (iday > 29) then
                    iday   = 1
                    imonth = imonth+1
                end if
            else
                if (iday > ndaycal(imonth,1)) then
                    iday   = 1
                    imonth = imonth+1
                end if
            end if

            if (imonth > 12) then
                imonth = 1
                iyear  = iyear+1
            end if
        
            ! additional variables to define forcing terms and boundary cond.
            if (iseasc >= 1) then
                imont1 = imonth
                tmonth = (iday-0.5)/float(ndaycal(imonth,1))
                tyear  = (ndaycal(imonth,2)+iday-0.5)/float(ncal)
            else
                imont1 = imonth0
                tmonth = 0.5
                tyear  = (ndaycal(imont1,2)&
                    & +0.5*ndaycal(imont1,2))/float(ncal)
            end if
        end subroutine newdate
end module
