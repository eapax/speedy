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
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: j1, j2
    real(dp), intent(in) :: dt, alph, rob, wil
    complex(dp), dimension(mx,nx,kx) :: vordt, divdt, tdt
    complex(dp) :: psdt(mx,nx), trdt(mx,nx,kx,ntr)
    real(dp) :: eps, sdrag

    complex(dp) :: ctmp(mx,nx,kx)

    integer :: iitest = 0, itr, k

    if (iitest==1) print*, ' inside step'

    ! 1. Computation of grid-point tendencies
    ! (converted to spectral at the end of GRTEND)
    if (iitest==1) print*,' call grtend'
    call grtend(vordt,divdt,tdt,psdt,trdt,1,j2)

    ! 2. Computation of spectral tendencies
    if (alph==0.0_dp) then
        if (iitest==1) print*,' call sptend'
        call sptend(divdt,tdt,psdt,j2)
    else
        if (iitest==1) print*,' call sptend'
        call sptend(divdt,tdt,psdt,1)

        ! implicit correction
        if (iitest==1) print*,' call implic'
        call implic(divdt,tdt,psdt)
    endif

    ! 3. Horizontal diffusion
    if (iitest==1) print*, ' biharmonic damping '

    ! 3.1 Diffusion of wind and temperature
    call hordif(kx,vor,vordt,dmp, dmp1)
    call hordif(kx,div,divdt,dmpd,dmp1d)

    do k=1,kx
        ctmp(:,:,k) = t(:,:,k,1) + tcorh*tcorv(k)
    enddo

    call hordif(kx,ctmp,tdt,dmp,dmp1)

    ! 3.2 Stratospheric diffusion and zonal wind damping
    sdrag = 1.0_dp/(tdrs*3600.0_dp)
    vordt(1,:,1) = vordt(1,:,1)-sdrag*vor(1,:,1,1)
    divdt(1,:,1) = divdt(1,:,1)-sdrag*div(1,:,1,1)

    call hordif(1,vor, vordt,dmps,dmp1s)
    call hordif(1,div, divdt,dmps,dmp1s)
    call hordif(1,ctmp,tdt,  dmps,dmp1s)

    ! 3.4 Diffusion of tracers
    do k=1,kx
        ctmp(:,:,k) = tr(:,:,k,1,1) + qcorh*qcorv(k)
    enddo

    call hordif(kx,ctmp,trdt,dmpd,dmp1d)

    if (ntr>1) then
        do itr=2,ntr
            call hordif(kx,tr(:,:,:,1,itr),trdt(:,:,:,itr),dmp,dmp1)
        enddo
    endif

    ! 4. Time integration with Robert filter
    if (dt<=0.0_dp) return

    if (iitest==1) print*,' time integration'
    if (j1==1) then
        eps = 0.0_dp
    else
        eps = rob
    endif

    call timint(j1,dt,eps,wil,1,ps,psdt)
    call timint(j1,dt,eps,wil,kx,vor,vordt)
    call timint(j1,dt,eps,wil,kx,div,divdt)
    call timint(j1,dt,eps,wil,kx,t,  tdt)

    do itr=1,ntr
        call timint(j1,dt,eps,wil,kx,tr(:,:,:,1,itr),trdt(:,:,:,itr))
    enddo
end subroutine step

subroutine hordif(nlev,field,fdt,dmp,dmp1)
    !   Aux. subr. HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
    !   Purpose : Add horizontal diffusion tendency of FIELD
    !             to spectral tendency FDT at NLEV levels
    !             using damping coefficients DMP and DMP1

    USE mod_atparam
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: nlev
    complex(dp), intent(in) :: field(mx,nx,kx)
    complex(dp), intent(inout) :: fdt(mx,nx,kx)
    real(dp), intent(in) :: dmp(mx,nx), dmp1(mx,nx)
    integer :: k, n, m

    do k=1,nlev
        do n=1,nx
            do m=1,mx
                fdt(m,n,k)=(fdt(m,n,k)-dmp(m,n)*field(m,n,k))*dmp1(m,n)
            end do
        end do
    end do
end subroutine hordif

subroutine timint(j1,dt,eps,wil,nlev,field,fdt)
    !  Aux. subr. timint (j1,dt,eps,wil,nlev,field,fdt)
    !  Purpose : Perform time integration of field at nlev levels
    !            using tendency fdt

    use mod_atparam
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: j1, nlev
    real(dp), intent(in) :: dt, eps, wil
    complex(dp), intent(in) :: fdt(mx,nx,nlev)
    complex(dp), intent(inout) :: field(mx,nx,nlev,2)
    real(dp) :: eps2
    complex(dp) :: fnew(mx,nx)
    integer :: k, n, m

    eps2 = 1.0_dp-2.0_dp*eps

    if (ix==iy*4) then
        do k=1,nlev
            call trunct(fdt(:,:,k))
        enddo
    endif

    ! The actual leap frog with the robert filter
    do k=1,nlev
        do n=1,nx
            do m=1,mx
                fnew (m,n)     = field(m,n,k,1) + dt*fdt(m,n,k)
                field(m,n,k,1) = field(m,n,k,j1) + wil*eps*(field(m,n,k,1)&
                    &-2.0_dp*field(m,n,k,j1)+fnew(m,n))

                ! and here comes Williams' innovation to the filter
                field(m,n,k,2) = fnew(m,n)-(1-wil)*eps*(field(m,n,k,1)&
                    &-2.0_dp*field(m,n,k,j1)+fnew(m,n))

            end do
        end do
    end do
end subroutine timint
