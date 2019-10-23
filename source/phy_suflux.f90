module phy_suflux
    use mod_atparam
    use humidity, only: shtorh_celsius, q_sat_celsius, zero_c
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_surface_fluxes, ini_suflux, truncate_suflux, suflux

    ! true : use an asymmetric stability coefficient
    logical, parameter :: lscasym = .true.
    ! true : use stability coef. to compute drag over sea
    logical, parameter :: lscdrag = .true.
    ! true : redefine skin temp. from energy balance
    logical, parameter :: lskineb = .true.

    ! Index for surface
    integer, parameter :: ks=2

    namelist /surface_fluxes/ fwind0, ftemp0, fhum0, cdl, cds, chl, chs, &
            vgust, ctday, dtheta, fstab, hdrag, fhdrag, clambda, clambsn

    !  Constants for surface fluxes
    ! Ratio of near-sfc wind to lowest-level wind
    type(rpe_var) :: fwind0

    ! Weight for near-sfc temperature extrapolation (0-1) :
    !          1 : linear extrapolation from two lowest levels
    !          0 : constant potential temperature ( = lowest level)
    type(rpe_var) :: ftemp0

    ! Weight for near-sfc specific humidity extrapolation (0-1) :
    !            1 : extrap. with constant relative hum. ( = lowest level)
    !            0 : constant specific hum. ( = lowest level)
    type(rpe_var) :: fhum0

    ! Drag coefficient for momentum over land
    type(rpe_var) :: cdl

    ! Drag coefficient for momentum over sea
    type(rpe_var) :: cds

    ! Heat exchange coefficient over land
    type(rpe_var) :: chl

    ! Heat exchange coefficient over sea
    type(rpe_var) :: chs

    ! Wind speed for sub-grid-scale gusts
    type(rpe_var) :: vgust

    ! Daily-cycle correction (dTskin/dSSRad)
    type(rpe_var) :: ctday

    ! Potential temp. gradient for stability correction
    type(rpe_var) :: dtheta

    ! Amplitude of stability correction (fraction)
    type(rpe_var) :: fstab

    ! Height scale for orographic correction
    type(rpe_var) :: hdrag

    ! Amplitude of orographic correction (fraction)
    type(rpe_var) :: fhdrag

    ! Heat conductivity in skin-to-root soil layer
    type(rpe_var) :: clambda

    ! Heat conductivity in soil for snow cover = 1
    type(rpe_var) :: clambsn

    ! Time-invariant fields (initial. in SFLSET)
    type(rpe_var), allocatable :: forog(:)

    ! Local derived variables
    type(rpe_var) :: esbc_1_3, esbc_1_4, prd, rcp, chlcp, chscp, rdphi0
    type(rpe_var), allocatable :: sqclat(:)

    ! Local copies of mod_physcon variables
    type(rpe_var) :: cp_sflx, alhc_sflx
    type(rpe_var), allocatable :: wvi_sflx(:,:)

    contains
        subroutine setup_surface_fluxes(fid)
            integer, intent(in) :: fid

            read(fid, surface_fluxes)

            allocate(forog(ngp))
            allocate(sqclat(il))

            write(*, surface_fluxes)

            allocate(wvi_sflx(kx,2))
        end subroutine setup_surface_fluxes

        subroutine ini_suflux()
            use mod_physcon, only: p0, rd, cp, alhc, sbc, clat, sigl, wvi
            use phy_radlw, only: emisfc

            ! Local derived variables
            esbc_1_3  = (emisfc*sbc)**(1.0_dp/3.0_dp)
            esbc_1_4  = (emisfc*sbc)**(1.0_dp/4.0_dp)
            prd = p0/rd
            rcp = 1.0_dp/cp
            chlcp = chl*cp
            chscp = chs*cp
            rdphi0 =-1.0_dp/(rd*288.0_dp*sigl(kx))
            sqclat=sqrt(clat(:))
            call sflset()

            ! Local copies of mod_physcon
            cp_sflx = cp
            alhc_sflx = alhc
            wvi_sflx = wvi
        end subroutine ini_suflux

        subroutine truncate_suflux()
            ! Namelist variables
            call apply_truncation(fwind0)
            call apply_truncation(ftemp0)
            call apply_truncation(fhum0)
            call apply_truncation(cdl)
            call apply_truncation(cds)
            call apply_truncation(chl)
            call apply_truncation(chs)
            call apply_truncation(vgust)
            call apply_truncation(ctday)
            call apply_truncation(dtheta)
            call apply_truncation(fstab)
            call apply_truncation(hdrag)
            call apply_truncation(fhdrag)
            call apply_truncation(clambda)
            call apply_truncation(clambsn)

            ! Local derived variables
            call apply_truncation(forog)
            call apply_truncation(esbc_1_3)
            call apply_truncation(esbc_1_4)
            call apply_truncation(prd)
            call apply_truncation(rcp)
            call apply_truncation(chlcp)
            call apply_truncation(chscp)
            call apply_truncation(rdphi0)
            call apply_truncation(sqclat)

            ! Local copies of mod_physcon
            call apply_truncation(cp_sflx)
            call apply_truncation(alhc_sflx)
            call apply_truncation(wvi_sflx)
        end subroutine truncate_suflux

        subroutine sflset()
            ! subroutine sflset ()
            !
            ! Purpose: compute orographic factor for land surface drag
            ! Input:   phi0   = surface geopotential            (2-dim)
            !          Initialized common blocks: sflfix

            use mod_surfcon, only: phis0
            use mod_physcon, only: gg

            integer :: i, j, ij
            type(rpe_var) :: rhdrag

            rhdrag = rpe_literal(1.0_dp)/(gg*hdrag)

            ij = 0
            do j = 1, il
                do i = 1, ix
                    ij = ij + 1
                    forog(ij) = rpe_literal(1.0_dp) + fhdrag*&
                            (rpe_literal(1.0_dp) - &
                                    exp(-max(phis0(i, j), rpe_literal(0.0_dp))*&
                                            rhdrag))

                end do
            end do
        end subroutine sflset

        subroutine suflux( &
                psa_in, ua_in, va_in, ta_in, qa_in, rh_in, phi_in, &
                phi0_in, fmask_in, tland_in, tsea_lm_in, tsea_om_in, &
                swav_in, ssrd_in, slrd_in, &
                hflx2tend_in, flx2tend_in, &
                ustr, vstr, shf, evap, slru, hfluxn, &
                tsfc, tskin, u0, v0, t0, q0, &
                ut_sflx, vt_sflx, tt_sflx, qt_sflx, &
                lfluxsea)
            !  Purpose: Compute surface fluxes of momentum, energy and moisture,
            !           and define surface skin temperature from energy balance

            ! The following variables from mod_fordate are calculated daily and
            ! then truncated to the precision for suflux as they are only used
            ! here.
            use mod_fordate, only: alb_l, snowc

            !  Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
            type(rpe_var), intent(in) :: psa_in(ngp)
            !           UA     = u-wind                          (3-dim)
            type(rpe_var), intent(in) :: ua_in(ngp,kx)
            !           VA     = v-wind                          (3-dim)
            type(rpe_var), intent(in) :: va_in(ngp,kx)
            !           TA     = temperature                     (3-dim)
            type(rpe_var), intent(in) :: ta_in(ngp,kx)
            !           QA     = specific humidity [g/kg]        (3-dim)
            type(rpe_var), intent(in) :: qa_in(ngp,kx)
            !           RH     = relative humidity [0-1]         (3-dim)
            type(rpe_var), intent(in) :: rh_in(ngp,kx)
            !           PHI    = geopotential                    (3-dim)
            type(rpe_var), intent(in) :: phi_in(ngp,kx)

            !           PHI0   = surface geopotential            (2-dim)
            type(rpe_var), intent(in) :: phi0_in(ngp)
            !           FMASK  = fractional land-sea mask        (2-dim)
            type(rpe_var), intent(in) :: fmask_in(ngp)
            !           TLAND  = land-surface temperature        (2-dim)
            type(rpe_var), intent(in) :: tland_in(ngp)
            !           TSEA_LM=  sea-surface temperature (land model) (2-dim)
            type(rpe_var), intent(in) :: tsea_lm_in(ngp)
            !           TSEA_OM=  sea-surface temperature (ocean model)(2-dim)
            type(rpe_var), intent(in) :: tsea_om_in(ngp)
            !           SWAV   = soil wetness availability [0-1] (2-dim)
            type(rpe_var), intent(in) :: swav_in(ngp)
            !           SSRD   = sfc sw radiation (downw. flux)  (2-dim)
            type(rpe_var), intent(in) :: ssrd_in(ngp)
            !           SLRD   = sfc lw radiation (downw. flux)  (2-dim)
            type(rpe_var), intent(in) :: slrd_in(ngp)

            ! hflx2tend = Conversion factor between heat fluxes and T tendency
            type(rpe_var), intent(in) :: hflx2tend_in(ngp,kx)
            ! flx2tend = Conversion factor between fluxes and tendencies
            type(rpe_var), intent(in) :: flx2tend_in(ngp,kx)

            !           LFLUXSEA   = Logical related to flux-correction
            logical, intent(in) :: lfluxsea

            !  Output:  USTR   = u stress                        (2-dim)
            type(rpe_var), intent(out) :: ustr(ngp,3)
            !           VSTR   = v stress                        (2-dim)
            type(rpe_var), intent(out) :: vstr(ngp,3)
            !           SHF    = sensible heat flux              (2-dim)
            type(rpe_var), intent(out) :: shf(ngp,3)
            !           EVAP   = evaporation [g/(m^2 s)]         (2-dim)
            type(rpe_var), intent(out) :: evap(ngp,3)
            !           SLRU   = sfc lw radiation (upward flux)  (2-dim)
            type(rpe_var), intent(out) :: slru(ngp,3)
            !           HFLUXN = net heat flux into land/sea     (2-dim)
            type(rpe_var), intent(out) :: hfluxn(ngp,2)

            !           TSFC   = surface temperature (clim.)     (2-dim)
            type(rpe_var), intent(out) :: tsfc(ngp)
            !           TSKIN  = skin surface temperature        (2-dim)
            type(rpe_var), intent(out) :: tskin(ngp)
            !           U0     = near-surface u-wind             (2-dim)
            type(rpe_var), intent(out) :: u0(ngp)
            !           V0     = near-surface v-wind             (2-dim)
            type(rpe_var), intent(out) :: v0(ngp)
            !           T0     = near-surface air temperature    (2-dim)
            type(rpe_var), intent(out) :: t0(ngp)
            !           Q0     = near-surface sp. humidity [g/kg](2-dim)
            type(rpe_var), intent(out) :: q0(ngp)

            ! ut_sflx = Zonal wind tendency due to surface fluxes
            type(rpe_var), intent(out) :: ut_sflx(ngp,kx)
            ! vt_sflx = Meridional wind tendency due to surface fluxes
            type(rpe_var), intent(out) :: vt_sflx(ngp,kx)
            ! tt_sflx = Temperature tendency due to surface fluxes
            type(rpe_var), intent(out) :: tt_sflx(ngp,kx)
            ! qt_sflx = Specific humidity tendency due to surface fluxes
            type(rpe_var), intent(out) :: qt_sflx(ngp,kx)

            ! Local copies of input variables
            type(rpe_var) :: psa(ngp)
            type(rpe_var), dimension(ngp,kx) :: ua, va, ta, qa, rh, phi
            type(rpe_var), dimension(ngp) :: &
                    phi0, fmask, tland, tsea_lm, tsea_om, swav, ssrd, slrd
            type(rpe_var), dimension(ngp,kx) :: hflx2tend, flx2tend

            ! Local variables
            integer :: j, j0, jlat

            type(rpe_var) :: denvvs(ngp, 0:2)
            type(rpe_var), dimension(ngp,2) :: t1, t2, q1, qsat0
            type(rpe_var) :: dslr(ngp), dtskin(ngp), clamb(ngp), astab, cdldv, &
                    cdsdv, dhfdt, dlambda, dt1, dthl, dths, &
                    ghum0, gtemp0, rdth, vg2

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            psa = psa_in
            ua = ua_in
            va = va_in
            ! Express Ta in Celsius. Should be done at double-precision then
            ! truncated as ta is not truncated on input
            ta = ta_in-zero_c
            qa = qa_in
            rh = rh_in
            phi = phi_in*rcp
            phi0 = phi0_in*rcp
            fmask = fmask_in
            ! Express other temperature quantities in Celsius
            tland = tland_in-zero_c
            tsea_lm = tsea_lm_in-zero_c
            tsea_om = tsea_om_in-zero_c
            swav = swav_in
            ssrd = ssrd_in
            slrd = slrd_in
            hflx2tend = hflx2tend_in
            flx2tend = flx2tend_in

            ! Initialisation
            ghum0 = rpe_literal(1.0_dp)-fhum0
            dlambda = clambsn-clambda

            ! 1. Extrapolation of wind, temp, hum. and density to the
            !    surface

            ! 1.1 Wind components
            u0 = fwind0 * ua(:,kx)
            v0 = fwind0 * va(:,kx)

            ! 1.2 Temperature
            gtemp0 = rpe_literal(1.0_dp)-ftemp0

            do j=1,ngp
                ! Temperature difference between lowest level and sfc
                dt1 = wvi_sflx(kx,2)*(ta(j,kx)-ta(j,kxm))
                ! Extrapolated temperature using actual lapse rate
                ! (1:land, 2:sea)
                t1(j,1) = ta(j,kx)+dt1
                t1(j,2) = t1(j,1)+phi0(j)*dt1*rdphi0*cp_sflx
                ! Extrapolated temperature using dry-adiab. lapse rate
                ! (1:land, 2:sea)
                t2(j,2) = ta(j,kx)+phi(j,kx)
                t2(j,1) = t2(j,2)-phi0(j)
            end do

            do j=1,ngp
                if (ta(j,kx)>ta(j,kxm)) then
                    ! Use extrapolated temp. if dT/dz < 0
                    t1(j,1) = ftemp0*t1(j,1)+gtemp0*t2(j,1)
                    t1(j,2) = ftemp0*t1(j,2)+gtemp0*t2(j,2)
                else
                    ! Use temp. at lowest level if dT/dz > 0
                    t1(j,1) = ta(j,kx)
                    t1(j,2) = ta(j,kx)
                endif
                t0(j) = t1(j,2)+fmask(j)*(t1(j,1)-t1(j,2))
            end do

            ! 1.3 Density * wind speed (including gustiness factor)
            vg2 = vgust*vgust

            do j=1,ngp
                denvvs(j,0)=(prd*psa(j)/(t0(j) + rpe_literal(zero_c)))*&
                        sqrt(u0(j)*u0(j) + v0(j)*v0(j) + vg2)
            end do

            ! 2. Compute land-sfc. fluxes using prescribed skin temperature

            ! 2.1 Define effective skin temperature to compensate for
            !     non-linearity of heat/moisture fluxes during the daily cycle
            do jlat=1,il
                j0=ix*(jlat-1)
                do j=j0+1,j0+ix
                    tskin(j)=tland(j) + ctday*sqclat(jlat)*ssrd(j)*&
                            (rpe_literal(1.0_dp)-alb_l(j))*psa(j)
                end do
            end do

            ! 2.2 Stability correction = f[pot.temp.(sfc)-pot.temp.(air)]
            rdth  = fstab/dtheta
            astab = 1.0_dp
            if (lscasym) astab = 0.5_dp
            ! to get smaller ds/dt in stable conditions

            do j=1,ngp
                ! Potential temp. difference (land+sea average)
                if (tskin(j)>t2(j,1)) then
                    dthl=min(dtheta,tskin(j)-t2(j,1))
                else
                    dthl=max(-dtheta,astab*(tskin(j)-t2(j,1)))
                endif
                denvvs(j,1)=denvvs(j,0)*(rpe_literal(1.0_dp)+dthl*rdth)
            end do

            ! 2.3 Wind stress
            do j=1,ngp
                cdldv     =  cdl*denvvs(j,0)*forog(j)
                ustr(j,1) = -cdldv*ua(j,kx)
                vstr(j,1) = -cdldv*va(j,kx)
            end do

            ! 2.4 Sensible heat flux
            do j=1,ngp
                shf(j,1) = chlcp*denvvs(j,1)*(tskin(j)-t1(j,1))
            end do

            ! 2.5 Evaporation
            if (fhum0>rpe_literal(0.0_dp)) then
                call shtorh_celsius(-1,ngp,t1(1,1),psa,rpe_literal(1.0_dp), &
                        q1(1,1),rh(1,kx),qsat0(1,1))

                do j=1,ngp
                  q1(j,1) = fhum0*q1(j,1)+ghum0*qa(j,kx)
                end do
            else
                q1(:,1) = qa(:,kx)
            end if

            qsat0(:, 1) = q_sat_celsius(ngp, tskin, psa, rpe_literal(1.0_dp))

            do j=1,ngp
                evap(j,1) = chl*denvvs(j,1) * &
                        max(rpe_literal(0.0_dp),swav(j)*qsat0(j,1)-q1(j,1))
            end do

            ! 3. Compute land-surface energy balance;
            !    adjust skin temperature and heat fluxes

            ! 3.1. Emission of lw radiation from the surface
            !      and net heat fluxes into land surface
            do j=1,ngp
                dslr(j)     = rpe_literal(4.0_dp)* &
                        (esbc_1_3*(tskin(j) + rpe_literal(zero_c)))**3
                slru(j,1)   = (esbc_1_4 *(tskin(j) + rpe_literal(zero_c)))**4
                hfluxn(j,1) = ssrd(j)*(rpe_literal(1.0_dp) - alb_l(j)) + &
                        slrd(j) - (slru(j,1) + shf(j,1) + alhc_sflx*evap(j,1))
            end do

            ! 3.2 Re-definition of skin temperature from energy balance
            if (lskineb) then
                ! Compute net heat flux including flux into ground
                do j=1,ngp
                  clamb(j)    = clambda+snowc(j)*dlambda
                  hfluxn(j,1) = hfluxn(j,1)-clamb(j)*(tskin(j)-tland(j))
                  dtskin(j)   = tskin(j)+rpe_literal(1.0_dp)
                end do

                ! Compute d(Evap) for a 1-degree increment of Tskin
                qsat0(:, 2) = q_sat_celsius(ngp, dtskin, psa, rpe_literal(1.0_dp))

                do j=1,ngp
                    if (evap(j,1)>0) then
                        qsat0(j,2) = swav(j)*(qsat0(j,2)-qsat0(j,1))
                    else
                        qsat0(j,2) = 0.0_dp
                    endif
                end do

                ! Redefine skin temperature to balance the heat budget
                do j=1,ngp
                    dhfdt = clamb(j) + dslr(j) + &
                            chl*denvvs(j,1)*(cp_sflx+alhc_sflx*qsat0(j,2))
                    dtskin(j) = hfluxn(j,1)/dhfdt
                    tskin(j)  = tskin(j)+dtskin(j)
                end do

                ! Add linear corrections to heat fluxes
                do j=1,ngp
                    shf(j,1)    = shf(j,1) +chlcp*denvvs(j,1)*dtskin(j)
                    evap(j,1)   = evap(j,1) + &
                            chl*denvvs(j,1)*qsat0(j,2)*dtskin(j)
                    slru(j,1)   = slru(j,1)+dslr(j)*dtskin(j)
                    hfluxn(j,1) = clamb(j)*(tskin(j)-tland(j))
                end do
            end if

            ! 4. Compute sea surface fluxes:
            ! 4.1 Correct near-sfc. air temperature over coastal sea points
            !     and compute near-sfc. humidity
            rdth  = fstab/dtheta
            astab = 1.0_dp
            if (lscasym) astab = 0.5_dp
            ! to get smaller dS/dT in stable conditions

            do j=1,ngp
                if (tsea_lm(j)>t2(j,2)) then
                   dths=min(dtheta,tsea_lm(j)-t2(j,2))
                else
                   dths=max(-dtheta,astab*(tsea_lm(j)-t2(j,2)))
                end if
                denvvs(j,2)=denvvs(j,0)*(rpe_literal(1.0_dp)+dths*rdth)
            end do

            if (fhum0>rpe_literal(0.0_dp)) then
                call shtorh_celsius(-1,ngp,t1(1,2),psa,rpe_literal(1.0_dp), &
                        q1(1,2),rh(1,kx),qsat0(1,2))

                do j=1,ngp
                  q1(j,2) = fhum0*q1(j,2)+ghum0*qa(j,kx)
                end do
            else
                q1(:,2) = qa(:,kx)
            end if

            ! 4.2 Wind stress
            do j=1,ngp
                cdsdv     =  cds*denvvs(j,ks)
                ustr(j,2) = -cdsdv*ua(j,kx)
                vstr(j,2) = -cdsdv*va(j,kx)
            end do

            ! Start of sea-sfc. heat fluxes computation
            call sea_fluxes(psa, tsea_lm, ssrd, slrd, &
                    denvvs(:,ks), t1(:,2), q1(:,2), &
                    shf(:,2), evap(:,2), slru(:,2), hfluxn(:,2))

            ! 3. Weighted average of surface fluxes and temperatures
            !    according to land-sea mask
            do j=1,ngp
                ustr(j,3) = ustr(j,2)+fmask(j)*(ustr(j,1)-ustr(j,2))
                vstr(j,3) = vstr(j,2)+fmask(j)*(vstr(j,1)-vstr(j,2))
                shf(j,3) =  shf(j,2)+fmask(j)*( shf(j,1)- shf(j,2))
                evap(j,3) = evap(j,2)+fmask(j)*(evap(j,1)-evap(j,2))
                slru(j,3) = slru(j,2)+fmask(j)*(slru(j,1)-slru(j,2))
            end do

            do j=1,ngp
                tsfc(j)  = tsea_lm(j)+fmask(j)*(tland(j)-tsea_lm(j))
                tskin(j) = tsea_lm(j)+fmask(j)*(tskin(j)-tsea_lm(j))
                t0(j)    = t1(j,2)+fmask(j)*(t1(j,1)- t1(j,2))
                q0(j)    = q1(j,2)+fmask(j)*(q1(j,1)- q1(j,2))
            end do
            ! End of 'land-mode' computation

            if (lfluxsea) then
                ! Recompute sea fluxes using sst from the ocean model in case of
                ! anomaly coupling
                call sea_fluxes(psa, tsea_om, ssrd, slrd, &
                        denvvs(:,ks), t1(:,2), q1(:,2), &
                        shf(:,2), evap(:,2), slru(:,2), hfluxn(:,2))
            end if

            ! Convert surface fluxes to tendencies
            do j=1,ngp
                ut_sflx(j,kx)=ustr(j,3)* flx2tend(j,kx)
                vt_sflx(j,kx)=vstr(j,3)* flx2tend(j,kx)
                tt_sflx(j,kx)= shf(j,3)*hflx2tend(j,kx)
                qt_sflx(j,kx)=evap(j,3)* flx2tend(j,kx)
            end do

            ! Convert output Temperatures back to Kelvin
            tsfc = tsfc + rpe_literal(zero_c)
            tskin = tskin + rpe_literal(zero_c)
        end subroutine suflux

        subroutine sea_fluxes(psa, tsea, ssrd, slrd, denvvs, t1, q1, &
                shf, evap, slru, hfluxn)
            ! Start of sea-sfc. heat fluxes computation

            ! The following variables from mod_fordate are calculated daily and
            ! then truncated to the precision for suflux as they are only used
            ! here.
            use mod_fordate, only: alb_s

            !  Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
            type(rpe_var), intent(in) :: psa(ngp)
            !           TSEA   =  sea-surface temperature        (2-dim)
            type(rpe_var), intent(in) :: tsea(ngp)
            !           SSRD   = sfc sw radiation (downw. flux)  (2-dim)
            type(rpe_var), intent(in) :: ssrd(ngp)
            !           SLRD   = sfc lw radiation (downw. flux)  (2-dim)
            type(rpe_var), intent(in) :: slrd(ngp)

            type(rpe_var), intent(in) :: denvvs(ngp)
            type(rpe_var), intent(in) :: t1(ngp)
            type(rpe_var), intent(in) :: q1(ngp)

            !  Output:  SHF    = sensible heat flux              (2-dim)
            type(rpe_var), intent(out) :: shf(ngp)
            !           EVAP   = evaporation [g/(m^2 s)]         (2-dim)
            type(rpe_var), intent(out) :: evap(ngp)
            !           SLRU   = sfc lw radiation (upward flux)  (2-dim)
            type(rpe_var), intent(out) :: slru(ngp)
            !           HFLUXN = net heat flux into land/sea     (2-dim)
            type(rpe_var), intent(out) :: hfluxn(ngp)

            ! Local variables
            type(rpe_var) :: qsat0(ngp)
            integer :: j

            ! 4.3 Sensible heat flux
            do j=1,ngp
                shf(j) = chscp*denvvs(j)*(tsea(j)-t1(j))
            end do

            ! 4.4 Evaporation
            qsat0 = q_sat_celsius(ngp, tsea, psa, rpe_literal(1.0_dp))

            do j=1,ngp
                evap(j) = chs*denvvs(j)*(qsat0(j)-q1(j))
            end do

            ! 4.5 Emission of lw radiation from the surface
            !     and net heat fluxes into sea surface
            do j=1,ngp
                slru(j)   = (esbc_1_4 *(tsea(j) + rpe_literal(zero_c)))**4
                hfluxn(j) = ssrd(j)*(rpe_literal(1.0_dp)-alb_s(j))+slrd(j)-&
                    & (slru(j)+shf(j)+alhc_sflx*evap(j))
            end do

        end subroutine sea_fluxes
end module phy_suflux
