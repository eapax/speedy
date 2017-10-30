subroutine stloop(istep)
    !   subroutine stloop (istep)
    !
    !   Purpose: Perform a series of time steps calling 
    !            post-processing/output routines at selected steps
    !   Input/output : istep = time step index
    !   Updated common block : lflag2
      
    use mod_lflags, only: lradsw, lrandf
    use mod_tsteps
    use mod_date, only: ihour, iday, update_time
    use mod_dynvar
    use ppo_IO_stream, only: IO_stream, init_IO_stream, write_IO_stream

    implicit none

    integer, intent(inout) :: istep
    integer :: iitest = 0, j

    type(IO_stream) :: test_stream

    test_stream = init_IO_stream(filename='spectra.grd', spectral=.true., &
            var_ID=(/ 1,2,3,5 /))
    ! Loop over number of steps per day
    do j = 1, nsteps
        if (iitest == 1) print*, 'stloop: calling step ', istep

        ! Increment the time by one timestep
        call update_time()

        ! Set logical flags
        lradsw = (mod(istep,nstrad) == 1)
        lrandf = ((istep <= nstrdf) .or. (nstrdf < 0))

        ! Perform one leapfrog time step
        call step(2, 2, delt2, alph, rob, wil)

        ! Do diagnostic, post-processing and I/O tasks
        call diagns(2, istep)

        call write_IO_stream(test_stream)

        if (sixhrrun .and. ihour.eq.6) then
            call restart (2)
            print *,'normal end with 6-hr fcst (yeahhhhhhh!!!!)'
            stop 1111
        end if
        istep = istep + 1
    end do

    stop

end
