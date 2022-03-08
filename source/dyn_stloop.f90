subroutine stloop(istep)
    !   subroutine stloop (istep)
    !
    !   Purpose: Perform a series of time steps calling
    !            post-processing/output routines at selected steps
    !   Input/output : istep = time step index
    !   Updated common block : lflag2


    use mod_tsteps, only: nsteps, delt2, alph, rob, wil
    use mod_date, only: update_time
    use phy_radsw, only: lradsw, nstrad
    use rp_emulator
    use mod_prec, only: set_precision


    implicit none

    integer, intent(inout) :: istep
    integer :: j

    ! Loop over number of steps per day
    do j = 1, nsteps
        ! Increment the time by one timestep
        call update_time()

        ! Set logical flags
        lradsw = (mod(istep,nstrad)==1)

        ! Perform one leapfrog time step
        !call set_precision('rp_step')
        call step(2, 2, delt2, alph, rob, wil)
        !call set_precision('rp_stloop')
        
        ! Do diagnostic, post-processing and I/O tasks
        call diagns(2, istep)

        istep = istep + 1

    end do
end subroutine
