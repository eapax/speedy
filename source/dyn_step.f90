subroutine step(j1,j2,dt,alph,rob,wil)
    !   subroutine step (j1,j2,dt,alph,rob,wil)
    !
    !   Purpose: perform one time step starting from F(1) and F(2) 
    !            and using the following scheme:
    !
    !   Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
    !   F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
    !   F(2) = Fnew
    !
    !   Input: 
    !   If J1=1, J2=1 : forward time step (eps=0)
    !   If J1=1, J2=2 : initial leapfrog time step (eps=0)
    !   If J1=2, J2=2 : leapfrog time step with time filter (eps=ROB)
    !   DT = time step (if DT < or = 0, tendencies are computed but 
    !                   no time stepping is performed)
    !   alph = 0   : forward step for gravity wave terms
    !   alph = 1   : backward implicit step for g.w.
    !   alph = 0.5 : centered implicit step for g.w.
    !   rob  = Robert filter coefficient
    !   wil  = Williams filter coefficient
     
    use mod_dyncon0, only: tdrs
    use mod_atparam
    use mod_dynvar
    use mod_hdifcon
    use rp_emulator
    use mod_prec, only: set_precision

    implicit none

    integer, intent(in) :: j1, j2
    type(rpe_var), intent(in) :: dt, alph, rob, wil
    type(rpe_complex_var), dimension(mx,nx,kx) :: vordt, divdt, tdt
    type(rpe_complex_var) :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
    type(rpe_var) :: eps, sdrag

    type(rpe_complex_var) :: ctmp(mx,nx,kx)

    integer :: iitest = 0, n, itr, k, m

    if (iitest.eq.1) print*, ' inside step'

    ! 1. Computation of grid-point tendencies
    ! (converted to spectral at the end of GRTEND)
    if (iitest.eq.1) print*,' call grtend'
    call grtend(vordt,divdt,tdt,psdt,trdt,1,j2)

    ! 2. Computation of spectral tendencies
    call set_precision('Spectral Dynamics')
    if (alph.eq.0.) then
        if (iitest.eq.1) print*,' call sptend'
        call sptend(divdt,tdt,psdt,j2)
    else
        if (iitest.eq.1) print*,' call sptend'
        call sptend(divdt,tdt,psdt,1)

        ! implicit correction 
        if (iitest.eq.1) print*,' call implic'
        call implic(divdt,tdt,psdt)
    endif

    ! 3. Horizontal diffusion
    call set_precision('Diffusion')
    if (iitest.eq.1) print*, ' biharmonic damping '

    ! 3.1 Diffusion of wind and temperature
    call hordif(kx,vor,vordt,dmp, dmp1)
    call hordif(kx,div,divdt,dmpd,dmp1d)

    do k=1,kx
        ctmp(:,:,k) = t(:,:,k,1) + tcorh*tcorv(k)
    enddo

    call hordif(kx,ctmp,tdt,dmp,dmp1)

    ! 3.2 Stratospheric diffusion and zonal wind damping
    sdrag = rpe_literal(1.)/(tdrs*rpe_literal(3600.))
    vordt(1,:,1) = vordt(1,:,1)-sdrag*vor(1,:,1,1)
    divdt(1,:,1) = divdt(1,:,1)-sdrag*div(1,:,1,1)

    call hordif(1,vor, vordt,dmps,dmp1s)
    call hordif(1,div, divdt,dmps,dmp1s)
    call hordif(1,ctmp,tdt,  dmps,dmp1s)

    ! 3.3 Check for eddy kinetic energy growth rate 
    ! CALL CGRATE (VOR,DIV,VORDT,DIVDT)

    ! 3.4 Diffusion of tracers
    do k=1,kx
        ctmp(:,:,k) = tr(:,:,k,1,1) + qcorh*qcorv(k)
    enddo

    call hordif(kx,ctmp,trdt,dmpd,dmp1d)

    if (ntr.gt.1) then
        do itr=2,ntr
            call hordif(kx,tr(1,1,1,1,itr),trdt(1,1,1,itr),dmp,dmp1)
        enddo
    endif

    ! 4. Time integration with Robert filter
    if (dt.le.0.) return

    if (iitest.eq.1) print*,' time integration'
    call set_precision('Tendencies')
    psdt = psdt
    vordt = vordt
    divdt = divdt
    tdt = tdt
    trdt = trdt

    call set_precision('Timestepping')

    if (j1.eq.1) then
        eps = 0.
    else
        eps = rob
    endif

    call timint(j1,dt,eps,wil,1,ps,psdt)
    call timint(j1,dt,eps,wil,kx,vor,vordt)
    call timint(j1,dt,eps,wil,kx,div,divdt)
    call timint(j1,dt,eps,wil,kx,t,  tdt)

    do itr=1,ntr
        call timint(j1,dt,eps,wil,kx,tr(1,1,1,1,itr),trdt(1,1,1,itr))
    enddo
end   

subroutine hordif(nlev,field,fdt,dmp,dmp1)
    !   Aux. subr. HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
    !   Purpose : Add horizontal diffusion tendency of FIELD 
    !             to spectral tendency FDT at NLEV levels
    !             using damping coefficients DMP and DMP1

    USE mod_atparam
    use rp_emulator

    implicit none

    integer, intent(in) :: nlev
    type(rpe_complex_var), intent(in) :: field(mxnx,kx)
    type(rpe_complex_var), intent(inout) :: fdt(mxnx,kx)
    type(rpe_var), intent(in) :: dmp(mxnx), dmp1(mxnx)
    integer :: k, m

    do k=1,nlev
        fdt(:,k)=(fdt(:,k)-dmp*field(:,k))*dmp1
    enddo
end

subroutine timint(j1,dt,eps,wil,nlev,field,fdt)
    !  Aux. subr. timint (j1,dt,eps,wil,nlev,field,fdt)
    !  Purpose : Perform time integration of field at nlev levels
    !            using tendency fdt

    use mod_atparam
    use rp_emulator

    implicit none

    integer, intent(in) :: j1, nlev
    type(rpe_var), intent(in) :: dt, eps, wil
    type(rpe_complex_var), intent(in) :: fdt(mxnx,nlev)
    type(rpe_complex_var), intent(inout) :: field(mxnx,nlev,2)
    type(rpe_var) :: eps2, two
    type(rpe_complex_var) :: fnew(mxnx)
    integer :: k, m

    two = 2.0

    eps2 = rpe_literal(1.)-two*eps

    if (ix.eq.iy*4) then
        do k=1,nlev
            call trunct(fdt(1,k))
        enddo
    endif

    ! The actual leap frog with the robert filter
    do k=1,nlev
        fnew = field(:,k,1) + dt*fdt(:,k)
        field(:,k,1) = field(:,k,j1) +  wil*eps*(field(:,k,1)&
                & -two*field(:,k,j1)+fnew)
        ! and here comes Williams' innovation to the filter
        field(:,k,2) = fnew - (1-wil)*eps* &
                (field(:,k,1) - two*field(:,k,j1)+fnew)
    enddo
end

subroutine cgrate(vor,div,vordt,divdt)
    !   SUBROUTINE CGRATE (VOR,DIV,VORDT,DIVDT)
    !
    !   Purpose: Check growth rate of eddy kin. energy 
    !   Input  : VOR    = vorticity
    !            DIV    = divergence
    !            VORDT  = time derivative of VOR
    !            DIVDT  = time derivative of DIV
    
    USE mod_atparam
    use spectral, only: invlap
    use rp_emulator

    implicit none

    type(rpe_complex_var), dimension(mx,nx,kx), intent(in) :: vor, div
    type(rpe_complex_var), dimension(mx,nx,kx), intent(inout) :: vordt, divdt
    type(rpe_complex_var) :: temp(mx,nx)
    type(rpe_var) :: cdamp, grate, grmax, rnorm
    integer :: k, m, n

    grmax=0.2/(86400.*2.)

    cdamp=0.

    do k=2,kx
        grate=0.
        rnorm=0.

        call invlap (vor(:,:,k),temp)

        do n=1,nx
            do m=2,mx
                grate=grate-realpart(vordt(m,n,k)*conjg(temp(m,n)))
                rnorm=rnorm-realpart(  vor(m,n,k)*conjg(temp(m,n)))
            enddo
        enddo

        if (grate.gt.grmax*rnorm) cdamp = max(cdamp,0.8*grate/rnorm)
        ! if (grate.gt.grmax*rnorm) cdamp =&
        !     & max(cdamp,(grate*grate)/(grmax*rnorm*rnorm))
    enddo

    if (cdamp.gt.0.) then
        print *, ' rot. wind damping enabled'

        do k=1,kx
            do n=1,nx
                do m=2,mx
                    vordt(m,n,k)=vordt(m,n,k)-cdamp*vor(m,n,k)
                enddo
            enddo
        enddo
    endif

    cdamp=0.

    do k=2,kx
        grate=0.
        rnorm=0.

        call invlap (div(:,:,k),temp)

        do n=1,nx
            do m=2,mx
                grate=grate-realpart(divdt(m,n,k)*conjg(temp(m,n)))
                rnorm=rnorm-realpart(  div(m,n,k)*conjg(temp(m,n)))
            enddo
        enddo

        if (grate.gt.grmax*rnorm) cdamp = max(cdamp,0.8*grate/rnorm)
        !if (grate.gt.grmax*rnorm) cdamp =&
        !    & max(cdamp,(grate*grate)/(grmax*rnorm*rnorm))
    enddo

    if (cdamp.gt.0.) then
      print *, ' div. wind damping enabled'

      do k=1,kx
        do n=1,nx
         do m=2,mx
           divdt(m,n,k)=divdt(m,n,k)-cdamp*div(m,n,k)
         enddo
        enddo
      enddo
    endif
end
