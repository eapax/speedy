subroutine grtend(vordt,divdt,tdt,psdt,trdt,j1,j2)
    !   subroutine grtend (vordt,divdt,tdt,psdt,trdt,j1,j2)
    !
    !   Purpose: compute non-linear tendencies in grid-point space
    !            from dynamics and physical parametrizations,
    !            and convert them to spectral tendencies
    !
    !   dF/dt = T_dyn(F(J2)) + T_phy(F(J1))
    !
    !   Input:  j1 = time level index for physical tendencies
    !           j2 = time level index for dynamical tendencies
    !   Output: vordt = spectral tendency of vorticity
    !           divdt = spectral tendency of divergence
    !           tdt   = spectral tendency of temperature
    !           psdt  = spectral tendency of log(p_s)
    !           trdt  = spectral tendency of tracers
    !   Externals: The gridpoint fields stored in mod_physvar are updated in
    !              this subroutine
    !

    use mod_atparam
    use mod_dynvar, only: vor, div, t, ps, tr, phi
    use mod_physvar, only: ug1, vg1, tg1, qg1, phig1, pslg1
    use mod_dyncon1, only: coriol
    use mod_physcon, only: cp
    use spectral, only: uvspec, grid
    use mod_prec, only: set_precision, dp
    use rp_emulator

    implicit none

    type(rpe_complex_var), dimension(mx,nx,kx),     intent(inout) :: vordt, &
            divdt, tdt
    type(rpe_complex_var), dimension(mx,nx),        intent(inout) :: psdt
    type(rpe_complex_var), dimension(mx,nx,kx,ntr), intent(inout) :: trdt
    integer, intent(in) :: j1, j2

    ! Gridpoint tendencies. Updated in phypar and dyntend
    type(rpe_var), dimension(ix,il,kx)     :: utend, vtend, ttend
    type(rpe_var), dimension(ix,il,kx,ntr) :: trtend

    ! Local gridpoint variables from spectral transform for calculating dynamics
    ! tendencies. Only calculate ug/vg/tg when j1/=j2
    type(rpe_var), dimension(ix,il,kx)     :: ug, vg, tg, vorg, divg
    type(rpe_var), dimension(ix,il,kx,ntr) :: trg

    integer :: k, i, itr, j
    !call set_precision('rp_gridfields')

    ! 1. Compute grid-point fields
    ! 1.1 Update geopotential in spectral space
    !call set_precision('rp_gridfields11')
    call geop(j1)
    !call set_precision('rp_gridfields')

    ! 1.2 Grid-point variables for physics tendencies
    !call set_precision('rp_gridfields12')
    do k=1,kx
      call uvspec(vor(:,:,k, j1),div(:,:,k, j1),ug1(:,k),vg1(:,k))
    end do

    ! Truncate variables where the spectral transform is still done at double
    ! precision
    ug1 = ug1 / 3600.0_dp
    vg1 = vg1 / 3600.0_dp

    do k=1,kx
      call grid(t(:,:,k,j1), tg1(:,k), 1)
    end do

    ! Normalise geopotential by cp to avoid overflows in physics
    do k=1,kx
      call grid(phi(:,:,k)*(1/cp), phig1(:,k), 1)
    end do

    ! Don't transform the two stratospheric levels where humidity is set to zero
    ! because it leads to overflows
    do k=3,kx
        call grid(tr(:,:,k,j1,1), qg1(:,k), 1)
    end do

    call grid(ps(:,:,j1),pslg1,1)
    !call set_precision('rp_gridfields')


    ! 1.3 Grid-point variables for dynamics tendencies
    !call set_precision('rp_gridfields13')

    ! Set units of vorticity and divergence to 'per hour' to reduce underflows
    do k=1,kx
        call grid(vor(:,:,k,j2),vorg(:,:,k),1)
        call grid(div(:,:,k,j2),divg(:,:,k),1)

        do itr=1,ntr
          call grid(tr(:,:,k,j2,itr),trg(:,:,k,itr),1)
        end do

        do j=1,il
            do i=1,ix
                vorg(i,j,k)=vorg(i,j,k)+coriol(j)
            end do
        end do
    end do
    !call set_precision('rp_gridfields')

    !call set_precision('rp_grtend')

    ! 2. Parametrized physics tendencies
    !call set_precision('rp_phypar')
    call phypar(utend, vtend, ttend, trtend)
    !call set_precision('rp_grtend')


    ! 3. Dynamics tendencies
    !call set_precision('rp_dyntend')

    if (j1==j2) then
        call dyntend(vordt, divdt, tdt, psdt, trdt, j2, &
                     utend, vtend, ttend, trtend, &
                     ug1, vg1, tg1, vorg, divg, trg)
    else
        ! Only recalculate tg, ug and vg if phypar and dyntend use different
        ! time indices (j1/=j2)
        do k=1,kx
            call uvspec(vor(:,:,k,j2), div(:,:,k,j2), ug(:,:,k), vg(:,:,k))
        end do

        ug = ug / 3600.0_dp
        vg = vg / 3600.0_dp

        do k=1,kx
            call grid(t(:,:,k,j2), tg(:,:,k),1)
        end do

        call dyntend(vordt, divdt, tdt, psdt, trdt, j2, &
                     utend, vtend, ttend, trtend, &
                     ug, vg, tg, vorg, divg, trg)
    end if


    !call set_precision('rp_grtend')



end subroutine grtend

subroutine dyntend(vordt, divdt, tdt, psdt, trdt, j2, &
                   utend, vtend, ttend, trtend, &
                   ug, vg, tg, vorg, divg, trg)
    !   Purpose: Compute non-linear tendencies in grid-point space from
    !            dynamics and add to physics tendencies. Convert total
    !            gridpoint tendencies to spectral tendencies
    !
    !   Input:  j2    = time level index for dynamical tendencies
    !           ug    = Gridpoint field of zonal velocity
    !           vg    = Gridpoint field of meridional velocity
    !           tg    = Gridpoint field of temperatur
    !           vorg  = Gridpoint field of vorticity
    !           divg  = Gridpoint field of divergence
    !           trg   = Gridpoint field of tracers
    !   InOut:  utend = gridpoint tendencity of zonal velocity
    !           vtend = gridpoint tendencity of meridional velocity
    !           ttend = gridpoint tendencity of temperature
    !           qtend = gridpoint tendencity of humidity
    !   Output: vordt = spectral tendency of vorticity
    !           divdt = spectral tendency of divergence
    !           tdt   = spectral tendency of temperature
    !           psdt  = spectral tendency of log(p_s)
    !           trdt  = spectral tendency of tracers
    use mod_atparam
    use mod_dynvar, only: ps
    use mod_dyncon1, only: akap, rgas, dhs, dhsr, fsgr
    use mod_dyncon2, only: tref, tref3
    use spectral, only: lap, grad, grid, spec, vdspec
    use humidity, only: zero_C
    use mod_prec, only: dp, set_precision
    use rp_emulator

    implicit none

    integer, intent(in) :: j2

    type(rpe_complex_var), dimension(mx,nx,kx),     intent(inout) :: vordt, &
            divdt, tdt
    type(rpe_complex_var), dimension(mx,nx),        intent(inout) :: psdt
    type(rpe_complex_var), dimension(mx,nx,kx,ntr), intent(inout) :: trdt

    type(rpe_var), dimension(ix,il,kx),     intent(inout) :: utend, vtend, ttend
    type(rpe_var), dimension(ix,il,kx,ntr), intent(inout) :: trtend

    type(rpe_var), dimension(ix,il,kx),     intent(in) :: ug, vg, tg, vorg, divg
    type(rpe_var), dimension(ix,il,kx,ntr), intent(in) :: trg

    type(rpe_complex_var) :: dumc(mx,nx,3)

    type(rpe_var), dimension(ix,il,kx) :: tgg, puv
    type(rpe_var), dimension(ix,il) :: px, py, umean, vmean, dmean
    type(rpe_var) :: sigdt(ix,il,kxp)
    type(rpe_var) :: temp(ix,il,kxp), sigm(ix,il,kxp), dumr(ix,il,3)

    integer :: i, j, k, itr

    umean(:,:) = 0.0_dp
    vmean(:,:) = 0.0_dp
    dmean(:,:) = 0.0_dp

    do k=1,kx
        umean = umean + ug(:,:,k) * dhs(k)
        vmean = vmean + vg(:,:,k) * dhs(k)
        dmean = dmean + divg(:,:,k) * dhs(k)
    end do

    ! Compute tendency of log(surface pressure)
    call grad(ps(:,:,j2),dumc(:,:,2),dumc(:,:,3))
    call grid(dumc(:,:,2)*3600.0_dp, px, 2)
    call grid(dumc(:,:,3)*3600.0_dp, py, 2)

    dumr(:,:,1) = -umean * px - vmean * py

    call spec(dumr(:,:,1),psdt)
    psdt(1,1)=(0.0_dp,0.0_dp)

    ! Compute "vertical" velocity
    sigdt(:,:,1) = 0.0_dp
    sigdt(:,:,kxp) = 0.0_dp
    sigm(:,:,1) = 0.0_dp
    sigm(:,:,kxp) = 0.0_dp

    ! (The following combination of terms is utilized later in the
    !     temperature equation)
    do k=1,kx
        puv(:,:,k) = (ug(:,:,k) - umean) * px + (vg(:,:,k) - vmean) * py
    end do

    do k=1,kx
        !cspj sigdt is the vertical velocity (in sigma coords)
        sigdt(:,:,k+1) = sigdt(:,:,k) - dhs(k)*(puv(:,:,k)+divg(:,:,k)-dmean)
        sigm(:,:,k+1) = sigm(:,:,k) - dhs(k)*puv(:,:,k)
    end do

    ! Subtract part of temperature field that is used as reference for
    ! implicit terms
    do k=1,kx
        tgg(:,:,k) = tg(:,:,k)-(tref(k)-zero_c)
    end do

    px = rgas*px
    py = rgas*py

    ! Zonal wind tendency
    temp(:,:,1) = 0.0_dp
    temp(:,:,kxp) = 0.0_dp

    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k) * (ug(:,:,k) - ug(:,:,k-1))
    end do

    do k=1,kx
        utend(:,:,k) = utend(:,:,k) + vg(:,:,k) * vorg(:,:,k) - &
                tgg(:,:,k)*px - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
    end do

    ! Meridional wind tendency
    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k) * (vg(:,:,k) - vg(:,:,k-1))
    end do

    do k=1,kx
        vtend(:,:,k) = vtend(:,:,k) - ug(:,:,k)*vorg(:,:,k) - &
                tgg(:,:,k)*py - (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
    end do

    ! Temperature tendency
    do k=2,kx
        temp(:,:,k) = sigdt(:,:,k)*(tgg(:,:,k) - tgg(:,:,k-1)) + &
                        sigm(:,:,k)*(tref(k) - tref(k-1))
    end do

    do k=1,kx
        do j=1,il
            do i=1,ix
                ttend(i,j,k)= ttend(i,j,k) + tgg(i,j,k)*divg(i,j,k)&
                    & -(temp(i,j,k+1)+temp(i,j,k))*dhsr(k)&
                    & +fsgr(k)*tgg(i,j,k)*(sigdt(i,j,k+1)+sigdt(i,j,k))&
                    & +tref3(k)*(sigm(i,j,k+1)+sigm(i,j,k))&
                    & +akap*((tg(i,j,k)+zero_c)*puv(i,j,k)&
                    & -tgg(i,j,k)*dmean(i,j))
            end do
        end do
    end do

    ! Tracer tendency

    do itr=1,ntr
        do k=2,kx
            temp(:,:,k)=sigdt(:,:,k)*(trg(:,:,k,itr)-trg(:,:,k-1,itr))
        end do

        !spj for moisture, vertical advection is not possible between top
        !spj two layers
        !kuch three layers
        !if(iinewtrace==1)then
        do k=2,3
            temp(:,:,k)=0.0_dp
        enddo
        !endif

        do k=1,kx
            trtend(:,:,k,itr)=trtend(:,:,k,itr) + trg(:,:,k,itr)*divg(:,:,k) - &
                    (temp(:,:,k+1) + temp(:,:,k))*dhsr(k)
        end do
    end do

    ! 4. Conversion of grid-point tendencies to spectral space and calculation
    !    of terms using grid-point and spectral components
    do k=1,kx
        !  convert u and v tendencies to vor and div spectral tendencies
        !  vdspec takes a grid u and a grid v and converts them to
        !  spectral vor and div

        call vdspec(utend(:,:,k),vtend(:,:,k),vordt(:,:,k),divdt(:,:,k),2)

        !  add lapl(0.5*(u**2+v**2)) to div tendency,
        !  and add div(vT) to spectral t tendency
        dumr(:,:,1)=rpe_literal(0.5_dp)* &
                (ug(:,:,k)*ug(:,:,k) + vg(:,:,k)*vg(:,:,k))*3600.0_dp
        dumr(:,:,2)=-ug(:,:,k)*tgg(:,:,k)*3600.0_dp
        dumr(:,:,3)=-vg(:,:,k)*tgg(:,:,k)*3600.0_dp

        !  divergence tendency
        call spec(dumr(:,:,1),dumc(:,:,1))
        call lap (dumc(:,:,1),dumc(:,:,2))

        !fk--   Change to keep dimensions
        divdt(:,:,k) = divdt(:,:,k) - dumc(:,:,2)

        !  temperature tendency
        call vdspec(dumr(:,:,2),dumr(:,:,3),dumc(:,:,1),tdt(:,:,k),2)
        call spec(ttend(:,:,k),dumc(:,:,2))

        !fk--   Change to keep dimensions
        tdt(:,:,k) = tdt(:,:,k) + dumc(:,:,2)

        ! tracer tendency
        do itr=1,ntr
            dumr(:,:,2)=-ug(:,:,k)*trg(:,:,k,itr)*3600.0_dp
            dumr(:,:,3)=-vg(:,:,k)*trg(:,:,k,itr)*3600.0_dp

            call spec(trtend(:,:,k,itr),dumc(:,:,2))
            call vdspec(dumr(:,:,2),dumr(:,:,3),dumc(:,:,1),trdt(:,:,k,itr),2)

            !fk--   Change to keep dimensions
            trdt(:,:,k,itr) = trdt(:,:,k,itr) + dumc(:,:,2)
        end do
    end do
end subroutine dyntend
