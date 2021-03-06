subroutine diagns(jj,istep)
    ! subroutine diagns(jj,istep)

    ! Purpose: print global means of eddy kinetic energy and temperature
    ! Input : jj    = time level index (1 or 2)
    !         istep = time step index


    use mod_tsteps, only: nstdia
    use mod_atparam
    use mod_dynvar, only: vor, div, t
    use spectral, only: invlap
    use ppo_output_stream, only: update_output
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: jj, istep

    integer :: k, m, n, kk
    type(rpe_complex_var) :: temp(mx,nx)
    real(dp) :: diag(kx,3), sqhalf

    call update_output(istep)

    ! 1. Get global-mean temperature and compute eddy kinetic energy
    sqhalf = sqrt(0.5_dp)

    do k=1,kx
        diag(k,1)=0.0_dp
        diag(k,2)=0.0_dp
        diag(k,3)=sqhalf*REAL(REAL(t(1,1,k,jj)%val))

        call invlap(vor(1,1,k,jj),temp)

        do m=2,mx
            do n=1,nx
                diag(k,1)=diag(k,1)-REAL(REAL(temp(m,n)%val*conjg(vor(m,n,k,jj)%val)))
            end do
        end do

        call invlap(div(1,1,k,jj),temp)

        do m=2,mx
            do n=1,nx
                diag(k,2)=diag(k,2)-REAL(REAL(temp(m,n)%val*conjg(div(m,n,k,jj)%val)))
            end do
        end do
    end do

    diag(:, 1) = diag(:, 1) * (1/3600.0_dp)**2
    diag(:, 2) = diag(:, 2) * (1/3600.0_dp)**2

    ! 2. Print results to screen
    if (mod(istep,nstdia)==0) then
        print 2001, istep, (diag(k,1),k=1,kx)
        print 2002,        (diag(k,2),k=1,kx)
        print 2003,        (diag(k,3),k=1,kx)
    end if

    ! 3. Stop integration if model variables are out of range
    do k=1,kx
        if (diag(k,1)>500 .or. diag(k,2)>500 .or. diag(k,3)<-100 .or.&
            & diag(k,3)>50) then

            print 2001, istep, (diag(kk,1),kk=1,kx)
            print 2002,        (diag(kk,2),kk=1,kx)
            print 2003,        (diag(kk,3),kk=1,kx)

            ! Dump model fields to restart file
            call restart(2)

            stop '*** model variables out of accepted range ***'
        end if
    end do

    2001 format(' step =',i6,' reke =',(10f8.2))
    2002 format         (13x,' deke =',(10f8.2))
    2003 format         (13x,' temp =',(10f8.2))
end subroutine diagns
