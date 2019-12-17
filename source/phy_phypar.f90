subroutine phypar(utend,vtend,ttend,qtend)
    ! Compute physical parametrization tendencies for u, v, t, q
    !
    ! Output arguments:
    !     utend: u-wind tendency (gp)
    !     vtend: v-wind tendency (gp)
    !     ttend: temp. tendency (gp)
    !     qtend: spec. hum. tendency (gp)
    !
    ! Modified common blocks:  mod_physvar

    use mod_atparam
    use mod_physvar
    use mod_cpl_flags, only: icsea
    use mod_physcon, only: sig, sigh, grdsig, grdscp, cp
    use mod_surfcon, only: fmask1, phis0
    use mod_var_land, only: stl_am, soilw_am
    use mod_var_sea, only: sst_am, ssti_om
    use mod_fluxes, only: increment_fluxes
    use humidity, only: shtorh_celsius, zero_c
    use phy_convmf, only: convmf
    use phy_lscond, only: lscond
    use phy_cloud, only: cloud
    use phy_radsw, only: lradsw, radsw
    use phy_radlw, only: radlw_transmissivity, radlw_down, radlw_up
    use phy_suflux, only: suflux
    use phy_vdifsc, only: vdifsc
    use phy_sppt, only: sppt_on, gen_sppt, additive_forcing
    use rp_emulator
    use mod_prec

    implicit none

    ! Physics tendencies of prognostic variables
    type(rpe_var), dimension(ngp,kx), intent(out) :: utend, vtend, ttend, qtend

    ! Index of the top layer of deep convection (diagnosed in convmf)
    integer :: iptop(ngp)

    ! Index of cloud top for radiation scheme (diagnosed in cloud)
    integer :: icltop(ngp), icnv(ngp)

    ! Conversion constant between (heat) fluxes and tendencies
    type(rpe_var) :: flx2tend(ngp,kx)
    type(rpe_var) :: hflx2tend(ngp,kx)

    integer :: j, k

    ! Reciprocal of surface pressure 1/psg
    type(rpe_var), dimension(ngp) :: rps

    ! gradient of dry static energy (dSE/dPHI)
    type(rpe_var), dimension(ngp) :: gse

    ! Temperature in Celsius
    tg1 = tg1 - zero_c

    ! Normalise geopotential by cp
!    phig1 = phig1 / cp

    ! Truncate all variables
    call set_precision('Half')
    call apply_truncation(pslg1)
    call apply_truncation(ug1)
    call apply_truncation(vg1)
    call apply_truncation(tg1)
    call apply_truncation(qg1)
    call apply_truncation(phig1)

    ! 1. Compute thermodynamic variables
    do j=1,ngp
        psg(j)=exp(pslg1(j))
        rps(j)=rpe_literal(1.0_dp)/psg(j)
    end do

    do k=1,kx
        do j=1,ngp
            ! Conversion for fluxes->tendencies
            flx2tend(j,k)  = rps(j)*grdsig(k)
            ! Conversion for heat fluxes->tendencies
            hflx2tend(j,k) = rps(j)*grdscp(k)
        end do
    end do

    do k=1,kx
        do j=1,ngp
            ! Remove negative humidity values
            qg1(j,k)=max(qg1(j,k),rpe_literal(0.0_dp))
        end do
    end do

    ! Use normalised geopotential and temperature in Celsius to calculate static
    ! energy as an anomaly
    se = tg1 + phig1

    do k=1,kx
        call shtorh_celsius(1,ngp,tg1(:,k),psg,sig(k),qg1(:,k),rh(:,k),qsat(:,k))
    end do

    ! 2. Precipitation
    ! 2.1 Deep convection
    call convmf(psg,se,qg1,qsat,flx2tend,&
            iptop,cbmf,precnv,tt_cnv,qt_cnv)

    do j=1,ngp
        icnv(j)=kx-iptop(j)
    end do

    ! 2.2 Large-scale condensation
    call lscond(psg,rh,qsat,iptop,precls,tt_lsc,qt_lsc)

    ! 3. Radiation (shortwave and longwave) and surface fluxes
    ! 3.1 Compute shortwave tendencies and initialize lw transmissivity
    ! The sw radiation may be called at selected time steps
    if (lradsw) then
        do j=1,ngp
            gse(j) = (se(j,kx-1)-se(j,kx))/(phig1(j,kx-1)-phig1(j,kx))
        end do

        call cloud(qg1,rh,precnv,precls,iptop,gse,fmask1,icltop,cloudc,clstr)

        do j=1,ngp
            cltop(j)=sigh(icltop(j)-1)*psg(j)
            prtop(j)=float(iptop(j))
        end do

        call radsw(psg,qg1,icltop,cloudc,clstr,hflx2tend,ssrd,ssr,tsr,tt_rsw)

        call radlw_transmissivity(psg, qg1, icltop, cloudc)
    end if

    ! 3.2 Compute downward longwave fluxes
    call radlw_down(tg1,slrd)

    ! 3.3. Compute surface fluxes and land skin temperature
    call suflux(psg,ug1,vg1,tg1,qg1,rh,phig1,&
            phis0,fmask1,stl_am,sst_am,ssti_om,soilw_am,ssrd,slrd,&
            hflx2tend, flx2tend, &
            ustr,vstr,shf,evap,slru,hfluxn,ts,tskin,u0,v0,t0,q0,&
            ut_sflx, vt_sflx, tt_sflx, qt_sflx, &
            icsea > 0)

    ! 3.4 Compute upward longwave fluxes, convert them to tendencies
    !     and add shortwave tendencies
    call radlw_up(tg1,ts,slrd,slru(:,3),hflx2tend,slr,olr,tt_rlw)

    ! 4. PBL interactions with lower troposphere
    ! 4.1 Vertical diffusion and shallow convection
    call vdifsc(se,rh,qg1,qsat,phig1,icnv,ut_pbl,vt_pbl,tt_pbl,qt_pbl)

    ! 5. Store all fluxes for coupling and daily-mean output
    call increment_fluxes()

    ! Sum physics tendencies
    ut_phy = ut_sflx
    vt_phy = vt_sflx
    tt_phy = tt_cnv + tt_lsc + tt_rsw + tt_sflx + tt_rlw + tt_pbl
    qt_phy = qt_cnv + qt_lsc + qt_sflx + qt_pbl

    call set_precision('Double')

    ! Add SPPT noise
    if (sppt_on) then
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

    ! Additive random noise
    call additive_forcing(tt_phy)

    ! Convert temperature back to Kelvin
    tg1 = tg1 + zero_c
end subroutine phypar
