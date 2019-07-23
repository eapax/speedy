module phy_suflux
    use mod_atparam
    use humidity, only: shtorh, q_sat
    use mod_prec, only: dp

    implicit none

    private
    public setup_surface_fluxes, ini_suflux, suflux

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
    real(dp) :: fwind0

    ! Weight for near-sfc temperature extrapolation (0-1) :
    !          1 : linear extrapolation from two lowest levels
    !          0 : constant potential temperature ( = lowest level)
    real(dp) :: ftemp0

    ! Weight for near-sfc specific humidity extrapolation (0-1) :
    !            1 : extrap. with constant relative hum. ( = lowest level)
    !            0 : constant specific hum. ( = lowest level)
    real(dp) :: fhum0

    ! Drag coefficient for momentum over land
    real(dp) :: cdl

    ! Drag coefficient for momentum over sea
    real(dp) :: cds

    ! Heat exchange coefficient over land
    real(dp) :: chl

    ! Heat exchange coefficient over sea
    real(dp) :: chs

    ! Wind speed for sub-grid-scale gusts
    real(dp) :: vgust

    ! Daily-cycle correction (dTskin/dSSRad)
    real(dp) :: ctday

    ! Potential temp. gradient for stability correction
    real(dp) :: dtheta

    ! Amplitude of stability correction (fraction)
    real(dp) :: fstab

    ! Height scale for orographic correction
    real(dp) :: hdrag

    ! Amplitude of orographic correction (fraction)
    real(dp) :: fhdrag

    ! Heat conductivity in skin-to-root soil layer
    real(dp) :: clambda

    ! Heat conductivity in soil for snow cover = 1
    real(dp) :: clambsn

    ! Time-invariant fields (initial. in SFLSET)
    real(dp), allocatable :: forog(:)

    ! Local derived variables
    real(dp) :: esbc, esbc4, prd, rcp, chlcp, chscp, rdphi0
    real(dp), allocatable :: sqclat(:)

    contains
        subroutine setup_surface_fluxes(fid)
            integer, intent(in) :: fid

            read(fid, surface_fluxes)

            allocate(forog(ngp))
            allocate(sqclat(il))

            write(*, surface_fluxes)
        end subroutine setup_surface_fluxes

        subroutine ini_suflux()
            use mod_physcon, only: p0, rd, cp, alhc, sbc, clat, sigl
            use phy_radlw, only: emisfc

            ! Local derived variables
            esbc  = emisfc*sbc
            esbc4 = 4.0_dp*esbc
            prd = p0/rd
            rcp = 1.0_dp/cp
            chlcp = chl*cp
            chscp = chs*cp
            rdphi0 =-1.0_dp/(rd*288.0_dp*sigl(kx))
            sqclat=sqrt(clat(:))
            call sflset()
        end subroutine ini_suflux

        subroutine sflset()
            ! subroutine sflset ()
            !
            ! Purpose: compute orographic factor for land surface drag
            ! Input:   phi0   = surface geopotential            (2-dim)
            !          Initialized common blocks: sflfix

            use mod_surfcon, only: phis0
            use mod_physcon, only: gg

            integer :: i, j, ij
            real(dp) :: rhdrag

            rhdrag = 1.0_dp/(gg*hdrag)

            ij = 0
            do j = 1, il
                do i = 1, ix
                    ij = ij + 1
                    forog(ij) = 1.0_dp + fhdrag*&
                            (1.0_dp - &
                                    exp(-max(phis0(i, j), 0.0_dp)*&
                                            rhdrag))

                end do
            end do
        end subroutine sflset

        subroutine suflux( &
                psa, ua, va, ta, qa, rh, phi, &
                phi0, fmask, tland, tsea_lm, tsea_om, swav, ssrd, slrd, &
                hflx2tend, flx2tend, &
                ustr, vstr, shf, evap, slru, hfluxn, &
                tsfc, tskin, u0, v0, t0, q0, &
                ut_sflx, vt_sflx, tt_sflx, qt_sflx, &
                lfluxsea)
            !  Purpose: Compute surface fluxes of momentum, energy and moisture,
            !           and define surface skin temperature from energy balance
            use mod_physcon, only: cp, alhc, wvi
            use mod_fordate, only: alb_l, snowc

            !  Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
            real(dp), intent(in) :: psa(ngp)
            !           UA     = u-wind                          (3-dim)
            real(dp), intent(in) :: ua(ngp,kx)
            !           VA     = v-wind                          (3-dim)
            real(dp), intent(in) :: va(ngp,kx)
            !           TA     = temperature                     (3-dim)
            real(dp), intent(in) :: ta(ngp,kx)
            !           QA     = specific humidity [g/kg]        (3-dim)
            real(dp), intent(in) :: qa(ngp,kx)
            !           RH     = relative humidity [0-1]         (3-dim)
            real(dp), intent(in) :: rh(ngp,kx)
            !           PHI    = geopotential                    (3-dim)
            real(dp), intent(in) :: phi(ngp,kx)

            !           PHI0   = surface geopotential            (2-dim)
            real(dp), intent(in) :: phi0(ngp)
            !           FMASK  = fractional land-sea mask        (2-dim)
            real(dp), intent(in) :: fmask(ngp)
            !           TLAND  = land-surface temperature        (2-dim)
            real(dp), intent(in) :: tland(ngp)
            !           TSEA_LM=  sea-surface temperature (land model) (2-dim)
            real(dp), intent(in) :: tsea_lm(ngp)
            !           TSEA_OM=  sea-surface temperature (ocean model)(2-dim)
            real(dp), intent(in) :: tsea_om(ngp)
            !           SWAV   = soil wetness availability [0-1] (2-dim)
            real(dp), intent(in) :: swav(ngp)
            !           SSRD   = sfc sw radiation (downw. flux)  (2-dim)
            real(dp), intent(in) :: ssrd(ngp)
            !           SLRD   = sfc lw radiation (downw. flux)  (2-dim)
            real(dp), intent(in) :: slrd(ngp)

            ! hflx2tend = Conversion factor between heat fluxes and T tendency
            real(dp), intent(in) :: hflx2tend(ngp,kx)
            ! flx2tend = Conversion factor between fluxes and tendencies
            real(dp), intent(in) :: flx2tend(ngp,kx)

            !           LFLUXSEA   = Logical related to flux-correction
            logical, intent(in) :: lfluxsea

            !  Output:  USTR   = u stress                        (2-dim)
            real(dp), intent(out) :: ustr(ngp,3)
            !           VSTR   = v stress                        (2-dim)
            real(dp), intent(out) :: vstr(ngp,3)
            !           SHF    = sensible heat flux              (2-dim)
            real(dp), intent(out) :: shf(ngp,3)
            !           EVAP   = evaporation [g/(m^2 s)]         (2-dim)
            real(dp), intent(out) :: evap(ngp,3)
            !           SLRU   = sfc lw radiation (upward flux)  (2-dim)
            real(dp), intent(out) :: slru(ngp,3)
            !           HFLUXN = net heat flux into land/sea     (2-dim)
            real(dp), intent(out) :: hfluxn(ngp,2)

            !           TSFC   = surface temperature (clim.)     (2-dim)
            real(dp), intent(out) :: tsfc(ngp)
            !           TSKIN  = skin surface temperature        (2-dim)
            real(dp), intent(out) :: tskin(ngp)
            !           U0     = near-surface u-wind             (2-dim)
            real(dp), intent(out) :: u0(ngp)
            !           V0     = near-surface v-wind             (2-dim)
            real(dp), intent(out) :: v0(ngp)
            !           T0     = near-surface air temperature    (2-dim)
            real(dp), intent(out) :: t0(ngp)
            !           Q0     = near-surface sp. humidity [g/kg](2-dim)
            real(dp), intent(out) :: q0(ngp)

            ! ut_sflx = Zonal wind tendency due to surface fluxes
            real(dp), intent(out) :: ut_sflx(ngp,kx)
            ! vt_sflx = Meridional wind tendency due to surface fluxes
            real(dp), intent(out) :: vt_sflx(ngp,kx)
            ! tt_sflx = Temperature tendency due to surface fluxes
            real(dp), intent(out) :: tt_sflx(ngp,kx)
            ! qt_sflx = Specific humidity tendency due to surface fluxes
            real(dp), intent(out) :: qt_sflx(ngp,kx)

            ! Local variables
            integer :: j, j0, jlat

            real(dp) :: denvvs(ngp, 0:2)
            real(dp), dimension(ngp,2) :: t1, t2, q1, qsat0
            real(dp) :: dslr(ngp), dtskin(ngp), clamb(ngp), astab, cdldv, &
                    cdsdv, dhfdt, dlambda, dt1, dthl, dths, &
                    ghum0, gtemp0, rdth, tsk3, vg2

            ! Initialisation
            ghum0 = 1.0_dp-fhum0
            dlambda = clambsn-clambda

            ! 1. Extrapolation of wind, temp, hum. and density to the
            !    surface

            ! 1.1 Wind components
            u0 = fwind0 * ua(:,kx)
            v0 = fwind0 * va(:,kx)

            ! 1.2 Temperature
            gtemp0 = 1.0_dp-ftemp0

            do j=1,ngp
                ! Temperature difference between lowest level and sfc
                dt1 = wvi(kx,2)*(ta(j,kx)-ta(j,kxm))
                ! Extrapolated temperature using actual lapse rate
                ! (1:land, 2:sea)
                t1(j,1) = ta(j,kx)+dt1
                t1(j,2) = t1(j,1)+phi0(j)*dt1*rdphi0
                ! Extrapolated temperature using dry-adiab. lapse rate
                ! (1:land, 2:sea)
                t2(j,2) = ta(j,kx)+rcp*phi(j,kx)
                t2(j,1) = t2(j,2)-rcp*phi0(j)
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
                denvvs(j,0)=(prd*psa(j)/t0(j))*sqrt(u0(j)*u0(j)+&
                        v0(j)*v0(j)+vg2)
            end do

            ! 2. Compute land-sfc. fluxes using prescribed skin temperature

            ! 2.1 Define effective skin temperature to compensate for
            !     non-linearity of heat/moisture fluxes during the daily cycle
            do jlat=1,il
                j0=ix*(jlat-1)
                do j=j0+1,j0+ix
                    tskin(j)=tland(j) + ctday*sqclat(jlat)*ssrd(j)*&
                            (1.0_dp-alb_l(j))*psa(j)
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
                denvvs(j,1)=denvvs(j,0)*(1.0_dp+dthl*rdth)
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
            if (fhum0>0.0_dp) then
                call shtorh(-1,ngp,t1(1,1),psa,1.0_dp, &
                        q1(1,1),rh(1,kx),qsat0(1,1))

                do j=1,ngp
                  q1(j,1) = fhum0*q1(j,1)+ghum0*qa(j,kx)
                end do
            else
                q1(:,1) = qa(:,kx)
            end if

            qsat0(:, 1) = q_sat(ngp, tskin, psa, 1.0_dp)

            do j=1,ngp
                evap(j,1) = chl*denvvs(j,1) * &
                        max(0.0_dp,swav(j)*qsat0(j,1)-q1(j,1))
            end do

            ! 3. Compute land-surface energy balance;
            !    adjust skin temperature and heat fluxes

            ! 3.1. Emission of lw radiation from the surface
            !      and net heat fluxes into land surface
            do j=1,ngp
                tsk3        = tskin(j)**3
                dslr(j)     = esbc4*tsk3
                slru(j,1)   = esbc *tsk3*tskin(j)
                hfluxn(j,1) = ssrd(j)*(1.0_dp - alb_l(j)) + &
                        slrd(j) - (slru(j,1) + shf(j,1) + alhc*evap(j,1))
            end do

            ! 3.2 Re-definition of skin temperature from energy balance
            if (lskineb) then
                ! Compute net heat flux including flux into ground
                do j=1,ngp
                  clamb(j)    = clambda+snowc(j)*dlambda
                  hfluxn(j,1) = hfluxn(j,1)-clamb(j)*(tskin(j)-tland(j))
                  dtskin(j)   = tskin(j)+1.0_dp
                end do

                ! Compute d(Evap) for a 1-degree increment of Tskin
                qsat0(:, 2) = q_sat(ngp, dtskin, psa, 1.0_dp)

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
                            chl*denvvs(j,1)*(cp+alhc*qsat0(j,2))
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
                denvvs(j,2)=denvvs(j,0)*(1.0_dp+dths*rdth)
            end do

            if (fhum0>0.0_dp) then
                call shtorh(-1,ngp,t1(1,2),psa,1.0_dp, &
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
        end subroutine suflux

        subroutine sea_fluxes(psa, tsea, ssrd, slrd, denvvs, t1, q1, &
                shf, evap, slru, hfluxn)
            ! Start of sea-sfc. heat fluxes computation

            ! The following variables from mod_fordate are calculated daily and
            ! only used here.
            use mod_physcon, only: alhc
            use mod_fordate, only: alb_s

            !  Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
            real(dp), intent(in) :: psa(ngp)
            !           TSEA   =  sea-surface temperature        (2-dim)
            real(dp), intent(in) :: tsea(ngp)
            !           SSRD   = sfc sw radiation (downw. flux)  (2-dim)
            real(dp), intent(in) :: ssrd(ngp)
            !           SLRD   = sfc lw radiation (downw. flux)  (2-dim)
            real(dp), intent(in) :: slrd(ngp)

            real(dp), intent(in) :: denvvs(ngp)
            real(dp), intent(in) :: t1(ngp)
            real(dp), intent(in) :: q1(ngp)

            !  Output:  SHF    = sensible heat flux              (2-dim)
            real(dp), intent(out) :: shf(ngp)
            !           EVAP   = evaporation [g/(m^2 s)]         (2-dim)
            real(dp), intent(out) :: evap(ngp)
            !           SLRU   = sfc lw radiation (upward flux)  (2-dim)
            real(dp), intent(out) :: slru(ngp)
            !           HFLUXN = net heat flux into land/sea     (2-dim)
            real(dp), intent(out) :: hfluxn(ngp)

            ! Local variables
            real(dp) :: qsat0(ngp)
            integer :: j

            ! 4.3 Sensible heat flux
            do j=1,ngp
                shf(j) = chscp*denvvs(j)*(tsea(j)-t1(j))
            end do

            ! 4.4 Evaporation
            qsat0 = q_sat(ngp, tsea, psa, 1.0_dp)

            do j=1,ngp
                evap(j) = chs*denvvs(j)*(qsat0(j)-q1(j))
            end do

            ! 4.5 Emission of lw radiation from the surface
            !     and net heat fluxes into sea surface
            do j=1,ngp
                slru(j)   = esbc*tsea(j)**4
                hfluxn(j) = ssrd(j)*(1.0_dp-alb_s(j))+slrd(j)-&
                    & (slru(j)+shf(j)+alhc*evap(j))
            end do

        end subroutine sea_fluxes
end module phy_suflux
