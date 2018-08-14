subroutine phypar(utend,vtend,ttend,qtend)
    !  subroutine phypar(utend,vtend,ttend,qtend)
    !
    !  Purpose: compute physical parametrization tendencies for u, v, t, q
    !  Output arguments:  utend  : u-wind tendency (gp)
    !                     vtend  : v-wind tendency (gp)
    !                     ttend  : temp. tendency (gp)
    !                     qtend  : spec. hum. tendency (gp)
    !  Modified common blocks:  mod_physvar

    use mod_atparam
    use mod_physvar
    use mod_cpl_flags, only: icsea
    use mod_physcon, only: sig, sigh, grdsig, grdscp, cp
    use mod_surfcon, only: fmask1, phis0
    use mod_var_land, only: stl_am, soilw_am
    use mod_var_sea, only: sst_am, ssti_om
    use mod_fluxes, only: increment_fluxes
    use humidity, only: shtorh
    use phy_convmf, only: convmf
    use phy_lscond, only: lscond
    use phy_cloud, only: cloud
    use phy_radsw, only: lradsw, radsw
    use phy_radlw, only: radlw
    use phy_suflux, only: suflux
    use phy_vdifsc, only: vdifsc
    use phy_sppt, only: sppt_on, gen_sppt
    use rp_emulator
    use mod_prec

    implicit none

    ! Physics tendencies of prognostic variables
    type(rpe_var), dimension(ngp,kx), intent(out) :: utend, vtend, ttend, qtend

    ! Index of the top layer of deep convection (diagnosed in convmf)
    integer :: iptop(ngp)

    ! Index of cloud top for radiation scheme (diagnosed in cloud)
    integer :: icltop(ngp), icnv(ngp)

    integer :: iitest=0, j, k

    ! Reciprocal of surface pressure 1/psg
    type(rpe_var), dimension(ngp) :: rps

    ! gradient of dry static energy (dSE/dPHI)
    type(rpe_var), dimension(ngp) :: gse

    ! 1. Compute thermodynamic variables
    if (iitest.eq.1) print *, ' 1.2 in phypar'

    do j=1,ngp
        psg(j)=exp(pslg1(j))
        rps(j)=rpe_literal(1.0_dp)/psg(j)
    end do

    do k=1,kx
        do j=1,ngp
            ! Remove negative humidity values
            qg1(j,k)=max(qg1(j,k),rpe_literal(0.0_dp))
            se(j,k)=cp*tg1(j,k)+phig1(j,k)
        end do
    end do

    do k=1,kx
        call shtorh(1,ngp,tg1(:,k),psg,sig(k),qg1(:,k),rh(:,k),qsat(:,k))
    end do

    ! 2. Precipitation
    ! 2.1 Deep convection
    call set_precision('Convection')
    call convmf(psg,se,qg1,qsat,iptop,cbmf,precnv,tt_cnv,qt_cnv)

    do k=2,kx
       do j=1,ngp
        tt_cnv(j,k) = tt_cnv(j,k)*rps(j)*grdscp(k)
        qt_cnv(j,k) = qt_cnv(j,k)*rps(j)*grdsig(k)
       end do
    end do

    do j=1,ngp
        icnv(j)=kx-iptop(j)
    end do

    ! 2.2 Large-scale condensation
    call set_precision('Condensation')
    call lscond(psg,qg1,qsat,iptop,precls,tt_lsc,qt_lsc)

    ! 3. Radiation (shortwave and longwave) and surface fluxes
    ! 3.1 Compute shortwave tendencies and initialize lw transmissivity
    if (iitest.eq.1) print *, ' 3.1 in PHYPAR'

    ! The sw radiation may be called at selected time steps
    if (lradsw) then
        call set_precision('Cloud')
        do j=1,ngp
            gse(j) = (se(j,kx-1)-se(j,kx))/(phig1(j,kx-1)-phig1(j,kx))
        end do

        call cloud(qg1,rh,precnv,precls,iptop,gse,fmask1,icltop,cloudc,clstr)

        do j=1,ngp
            cltop(j)=sigh(icltop(j)-1)*psg(j)
            prtop(j)=float(iptop(j))
        end do

        call set_precision('Short-Wave Radiation')
        call radsw(psg,qg1,icltop,cloudc,clstr,ssrd,ssr,tsr,tt_rsw)

        do k=1,kx
            do j=1,ngp
                tt_rsw(j,k)=tt_rsw(j,k)*rps(j)*grdscp(k)
            end do
        end do
    end if

    ! 3.2 Compute downward longwave fluxes
    call set_precision('Long-Wave Radiation')
    call radlw(-1,tg1,ts,slrd,slru(1,3),slr,olr,tt_rlw)

    ! 3.3. Compute surface fluxes and land skin temperature
    if (iitest.eq.1) then
        print *, ' 3.3 in PHYPAR'
        print *, 'mean(STL_AM) =', sum(STL_AM(:))/ngp
        print *, 'mean(SST_AM) =', sum(SST_AM(:))/ngp
    end if

    call set_precision('Surface Fluxes')
    call suflux(psg,ug1,vg1,tg1,qg1,rh,phig1,phis0,fmask1,stl_am,sst_am,&
        & soilw_am,ssrd,slrd,ustr,vstr,shf,evap,slru,hfluxn,ts,tskin,u0,v0,t0,&
        & q0,.true.)

    ! 3.3.1. Recompute sea fluxes in case of anomaly coupling
    if (icsea .gt. 0) then
       call suflux(psg,ug1,vg1,tg1,qg1,rh,phig1,phis0,fmask1,stl_am,ssti_om,&
           & soilw_am,ssrd,slrd,ustr,vstr,shf,evap,slru,hfluxn,ts,tskin,u0,v0,&
           & t0,q0,.false.)
    end if

    ! 3.4 Compute upward longwave fluxes, convert them to tendencies
    !     and add shortwave tendencies
    if (iitest.eq.1) print *, ' 3.4 in PHYPAR'
    call set_precision('Long-Wave Radiation')
    call radlw (1,tg1,ts,slrd,slru(1,3),slr,olr,tt_rlw)

    do k=1,kx
        do j=1,ngp
            tt_rlw(j,k) = tt_rlw(j,k)*rps(j)*grdscp(k)
        end do
    end do

    ! 4. PBL interactions with lower troposphere
    ! 4.1 Vertical diffusion and shallow convection
    call set_precision('Vertical Diffusion')
    call vdifsc(ug1,vg1,se,rh,qg1,qsat,phig1,icnv,ut_pbl,vt_pbl,tt_pbl,qt_pbl)

    ! 4.2 Add tendencies due to surface fluxes
    do j=1,ngp
        ut_pbl(j,kx)=ut_pbl(j,kx)+ustr(j,3)*rps(j)*grdsig(kx)
        vt_pbl(j,kx)=vt_pbl(j,kx)+vstr(j,3)*rps(j)*grdsig(kx)
        tt_pbl(j,kx)=tt_pbl(j,kx)+ shf(j,3)*rps(j)*grdscp(kx)
        qt_pbl(j,kx)=qt_pbl(j,kx)+evap(j,3)*rps(j)*grdsig(kx)
    end do

    ! 5. Store all fluxes for coupling and daily-mean output
    call set_precision('Grid Physics')
    call increment_fluxes()

    ! Sum physics tendencies
    ut_phy = ut_pbl
    vt_phy = vt_pbl
    tt_phy = tt_cnv + tt_lsc + tt_rsw + tt_rlw + tt_pbl
    qt_phy = qt_cnv + qt_lsc + qt_pbl

    ! Add SPPT noise
    if (sppt_on) then
        call set_precision('SPPT')
        sppt = gen_sppt()

        utend = sppt * ut_phy
        vtend = sppt * vt_phy
        ttend = sppt * tt_phy
        qtend = sppt * qt_phy

        ! Store SPPT tendencies
        ut_sppt = utend - ut_phy
        vt_sppt = vtend - vt_phy
        tt_sppt = ttend - tt_phy
        qt_sppt = qtend - qt_phy

    else
        utend = ut_phy
        vtend = vt_phy
        ttend = tt_phy
        qtend = qt_phy
    end if
end
