module surface_fluxes
    use mod_atparam
    use humidity, only: shtorh, q_sat
    use rp_emulator
    use mod_prec

    implicit none

    private
    public suflux, sflset, init_sflcon

    !  Constants for surface fluxes
    ! Ratio of near-sfc wind to lowest-level wind
    real(dp), parameter :: fwind0_ = 0.95_dp

    ! Weight for near-sfc temperature extrapolation (0-1) :
    !          1 : linear extrapolation from two lowest levels
    !          0 : constant potential temperature ( = lowest level)
    real(dp), parameter :: ftemp0_ = 1.0_dp

    ! Weight for near-sfc specific humidity extrapolation (0-1) :
    !            1 : extrap. with constant relative hum. ( = lowest level)
    !            0 : constant specific hum. ( = lowest level)
    real(dp), parameter :: fhum0_ = 0.0_dp

    ! Drag coefficient for momentum over land
    real(dp), parameter :: cdl_ = 2.4d-3

    ! Drag coefficient for momentum over sea
    real(dp), parameter :: cds_ = 1.0d-3

    ! Heat exchange coefficient over land
    real(dp), parameter :: chl_ = 1.2d-3

    ! Heat exchange coefficient over sea
    real(dp), parameter :: chs_ = 0.9d-3

    ! Wind speed for sub-grid-scale gusts
    real(dp), parameter :: vgust_ = 5.0_dp

    ! Daily-cycle correction (dTskin/dSSRad)
    real(dp), parameter :: ctday_ = 1.0d-2

    ! Potential temp. gradient for stability correction
    real(dp), parameter :: dtheta_ = 3.0_dp

    ! Amplitude of stability correction (fraction)
    real(dp), parameter :: fstab_ = 0.67_dp

    ! Height scale for orographic correction
    real(dp), parameter :: hdrag_ = 2000.0_dp

    ! Amplitude of orographic correction (fraction)
    real(dp), parameter :: fhdrag_ = 0.5_dp

    ! Heat conductivity in skin-to-root soil layer
    real(dp), parameter :: clambda_ = 7.0_dp

    ! Heat conductivity in soil for snow cover = 1
    real(dp), parameter :: clambsn_ = 7.0_dp

    ! Time-invariant fields (initial. in SFLSET)
    type(rpe_var) :: forog(ix*il)

    ! Reduced precision versions
    type(rpe_var) :: fwind0
    type(rpe_var) :: ftemp0
    type(rpe_var) :: fhum0
    type(rpe_var) :: cdl
    type(rpe_var) :: cds
    type(rpe_var) :: chl
    type(rpe_var) :: chs
    type(rpe_var) :: vgust
    type(rpe_var) :: ctday
    type(rpe_var) :: dtheta
    type(rpe_var) :: fstab
    type(rpe_var) :: hdrag
    type(rpe_var) :: fhdrag
    type(rpe_var) :: clambda
    type(rpe_var) :: clambsn

    contains

        subroutine init_sflcon
            fwind0  = fwind0_
            ftemp0  = ftemp0_
            fhum0   = fhum0_
            cdl     = cdl_
            cds     = cds_
            chl     = chl_
            chs     = chs_
            vgust   = vgust_
            ctday   = ctday_
            dtheta  = dtheta_
            fstab   = fstab_
            hdrag   = hdrag_
            fhdrag  = fhdrag_
            clambda = clambda_
            clambsn = clambsn_
        end subroutine

        subroutine suflux(psa,ua,va,ta,qa,rh,phi,phi0,fmask,tland,tsea,swav, &
                ssrd,slrd,ustr,vstr,shf,evap,slru,hfluxn,tsfc,tskin,u0,v0,t0, &
                q0,lfluxland)
            !  subroutine suflux (psa,ua,va,ta,qa,rh,phi,
            ! &                   phi0,fmask,tland,tsea,swav,ssrd,slrd,
            ! &                   ustr,vstr,shf,evap,slru,hfluxn,
            ! &                   tsfc,tskin,u0,v0,t0,q0,lfluxland)
            !
            !  Purpose: Compute surface fluxes of momentum, energy and moisture,
            !           and define surface skin temperature from energy balance
            !  Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
            !           UA     = u-wind                          (3-dim)
            !           VA     = v-wind                          (3-dim)
            !           TA     = temperature                     (3-dim)
            !           QA     = specific humidity [g/kg]        (3-dim)
            !           RH     = relative humidity [0-1]         (3-dim)
            !           PHI    = geopotential                    (3-dim)
            !           PHI0   = surface geopotential            (2-dim)
            !           FMASK  = fractional land-sea mask        (2-dim)
            !           TLAND  = land-surface temperature        (2-dim)
            !           TSEA   =  sea-surface temperature        (2-dim)
            !           SWAV   = soil wetness availability [0-1] (2-dim)
            !           SSRD   = sfc sw radiation (downw. flux)  (2-dim)
            !           SLRD   = sfc lw radiation (downw. flux)  (2-dim)
            !           LFLUXLAND   = Logical related ti flux-correction
            !  Output:  USTR   = u stress                        (2-dim)
            !           VSTR   = v stress                        (2-dim)
            !           SHF    = sensible heat flux              (2-dim)
            !           EVAP   = evaporation [g/(m^2 s)]         (2-dim)
            !           SLRU   = sfc lw radiation (upward flux)  (2-dim)
            !           HFLUXN = net heat flux into land/sea     (2-dim)
            !           TSFC   = surface temperature (clim.)     (2-dim)
            !           TSKIN  = skin surface temperature        (2-dim)
            !           U0     = near-surface u-wind             (2-dim)
            !           V0     = near-surface v-wind             (2-dim)
            !           T0     = near-surface air temperature    (2-dim)
            !           Q0     = near-surface sp. humidity [g/kg](2-dim)

            use mod_physcon, only: p0, rd, cp, alhc, sbc, sigl, wvi, clat
            use mod_radcon, only: emisfc, alb_l, alb_s, snowc
            use rp_emulator
            use mod_prec

            integer, parameter :: ngp=ix*il

            type(rpe_var), dimension(ngp,kx), intent(in) :: ua, va, ta, qa, rh, phi
            type(rpe_var), dimension(ngp), intent(in) :: phi0, fmask, tland, tsea, &
                    swav, ssrd, slrd

            type(rpe_var), dimension(ngp,3), intent(inout) :: ustr, vstr, shf, evap, slru
            type(rpe_var), intent(inout) :: hfluxn(ngp,2)
            type(rpe_var), dimension(ngp), intent(inout) :: tsfc, tskin, u0, v0, t0, q0

            integer :: j, j0, jlat, ks, nl1
            type(rpe_var), dimension(ngp,2), save :: t1, q1
            type(rpe_var), dimension(ngp,2) :: t2, qsat0
            type(rpe_var), save :: denvvs(ngp,0:2)
            type(rpe_var) :: dslr(ngp), dtskin(ngp), clamb(ngp), astab, cdldv, &
                    cdsdv, chlcp, chscp, dhfdt, dlambda, dt1, dthl, dths, &
                    esbc, esbc4, ghum0, gtemp0, prd, rcp, rdphi0, rdth, &
                    sqclat, tsk3, vg2

            logical lscasym, lscdrag, lskineb
            logical lfluxland

            type(rpe_var) :: psa(ngp)

            lscasym = .true.   ! true : use an asymmetric stability coefficient
            lscdrag = .true.   ! true : use stability coef. to compute drag
                               !        over sea
            lskineb = .true.   ! true : redefine skin temp. from energy balance

            !clambda = 7.       ! Heat conductivity in skin layer
            !clambsn = 7.       ! Heat conductivity for snow cover = 1

            esbc  = emisfc*sbc
            esbc4 = rpe_literal(4.)*esbc

            ghum0 = rpe_literal(1.)-fhum0

            dlambda = clambsn-clambda

            if (lfluxland)  then
                ! 1. Extrapolation of wind, temp, hum. and density to the
                !    surface

                ! 1.1 Wind components
                u0 = fwind0 * ua(:,kx)
                v0 = fwind0 * va(:,kx)

                ! 1.2 Temperature
                gtemp0 = rpe_literal(1.)-ftemp0
                rcp = rpe_literal(1.)/cp
                rdphi0 =-rpe_literal(1.)/(rd*rpe_literal(288.)*sigl(kx))
                nl1=kx-1

                do j=1,ngp
                    ! Temperature difference between lowest level and sfc
                    dt1 = wvi(kx,2)*(ta(j,kx)-ta(j,nl1))
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
                    if (ta(j,kx).gt.ta(j,nl1)) then
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

                ! 1.3 Spec. humidity
                !ghum0 = 1.-fhum0

                !call shtorh(-1,ngp,t0,psa,1.,q0,rh(1,kx),qsat0)

                !do j=1,ngp
                !    q0(j)=fhum0*q0(j)+ghum0*qa(j,kx)
                !end do

                ! 1.3 Density * wind speed (including gustiness factor)
                prd = p0/rd
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
                    sqclat=sqrt(clat(jlat))
                    do j=j0+1,j0+ix
                        tskin(j)=tland(j)+ &
                                ctday*sqclat*ssrd(j)*(1.-alb_l(j))*psa(j)
                    end do
                end do

                ! 2.2 Stability correction = f[pot.temp.(sfc)-pot.temp.(air)]
                rdth  = fstab/dtheta
                astab = 1.
                if (lscasym) astab = 0.5
                ! to get smaller ds/dt in stable conditions

                do j=1,ngp
                    ! Potential temp. difference (land+sea average)
                    !fkdth0 = tsea(j)-t2(j,2)
                    !fkdth0 = dth0+fmask(j)*((tskin(j)-t2(j,1))-dth0)

                    !fkif (dth0.gt.0.0) then
                    !fk   dthl=min(dtheta,dth0)
                    !fkelse
                    !fk   dthl=max(-dtheta,astab*dth0)
                    !fkendif

                    !fkdenvvs(j,1)=denvvs(j,0)*(1.+dthl*rdth)

                    if (tskin(j).gt.t2(j,1)) then
                        dthl=min(dtheta,tskin(j)-t2(j,1))
                    else
                        dthl=max(-dtheta,astab*(tskin(j)-t2(j,1)))
                    endif
                    denvvs(j,1)=denvvs(j,0)*(1.+dthl*rdth)
                end do

                ! 2.3 Wind stress
                do j=1,ngp
                    cdldv     =  cdl*denvvs(j,0)*forog(j)
                    ustr(j,1) = -cdldv*ua(j,kx)
                    vstr(j,1) = -cdldv*va(j,kx)
                end do

                ! 2.4 Sensible heat flux
                chlcp = chl*cp

                do j=1,ngp
                    shf(j,1) = chlcp*denvvs(j,1)*(tskin(j)-t1(j,1))
                end do

                ! 2.5 Evaporation
                if (fhum0.gt.0.) then
                    call shtorh(-1,ngp,t1(1,1),psa,rpe_literal(1.0), &
                            q1(1,1),rh(1,kx),qsat0(1,1))

                    do j=1,ngp
                      q1(j,1) = fhum0*q1(j,1)+ghum0*qa(j,kx)
                    end do
                else
                    q1(:,1) = qa(:,kx)
                end if

                qsat0(:, 1) = q_sat(ngp, tskin, psa, rpe_literal(1.0))

                do j=1,ngp
                    evap(j,1) = chl*denvvs(j,1) * &
                            max(0.0_dp,swav(j)*qsat0(j,1)-q1(j,1))
                end do

                ! 3. Compute land-surface energy balance;
                !    adjust skin temperature and heat fluxes

                ! 3.1. Emission of lw radiation from the surface
                !      and net heat fluxes into land surface
                do j=1,ngp
                    tsk3        = tskin(j)**rpe_literal(3)
                    dslr(j)     = esbc4*tsk3
                    slru(j,1)   = esbc *tsk3*tskin(j)
                    hfluxn(j,1) = ssrd(j)*(1.-alb_l(j))+slrd(j)-&
                        & (slru(j,1)+shf(j,1)+alhc*evap(j,1))
                end do

                ! 3.2 Re-definition of skin temperature from energy balance
                if (lskineb) then
                    ! Compute net heat flux including flux into ground
                    do j=1,ngp
                      clamb(j)    = clambda+snowc(j)*dlambda
                      hfluxn(j,1) = hfluxn(j,1)-clamb(j)*(tskin(j)-tland(j))
                      dtskin(j)   = tskin(j)+1.
                    end do

                    ! Compute d(Evap) for a 1-degree increment of Tskin
                    qsat0(:, 2) = q_sat(ngp, dtskin, psa, rpe_literal(1.0))

                    do j=1,ngp
                        if (evap(j,1).gt.0) then
                            qsat0(j,2) = swav(j)*(qsat0(j,2)-qsat0(j,1))
                        else
                            qsat0(j,2) = 0.
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
        !      ENDIF

                ! 4. Compute sea surface fluxes:
                !    Note: stability terms and wind stress are NOT re-defined
                !          if LFLUXLAND = .false.

                ! 4.1 Correct near-sfc. air temperature over coastal sea points
                !     and compute near-sfc. humidity
                !fkdo j=1,ngp
                !fk    if (fmask(j).gt.0.) then
                !fk        dtsea  = tsea(j) -t1(j,2)
                !fk        dtland = tskin(j)-t1(j,1)
                !fk        if (dtsea.gt.0.0.and.dtland.lt.0.0) then
                !fk            dtsea   = dtsea*(1.-fmask(j)**2)
                !fk            t1(j,2) = tsea(j)-dtsea
                !fk        endif
                !fk    endif
                !fkend do

                rdth  = fstab/dtheta
                astab = 1.
                if (lscasym) astab = 0.5
                ! to get smaller dS/dT in stable conditions

                do j=1,ngp
                    if (tsea(j).gt.t2(j,2)) then
                       dths=min(dtheta,tsea(j)-t2(j,2))
                    else
                       dths=max(-dtheta,astab*(tsea(j)-t2(j,2)))
                    end if
                    denvvs(j,2)=denvvs(j,0)*(1.+dths*rdth)
                end do

                if (fhum0.gt.0.) then
                    call shtorh(-1,ngp,t1(1,2),psa,rpe_literal(1.0), &
                            q1(1,2),rh(1,kx),qsat0(1,2))

                    do j=1,ngp
                      q1(j,2) = fhum0*q1(j,2)+ghum0*qa(j,kx)
                    end do
                else
                    q1(:,2) = qa(:,kx)
                end if

                ! 4.2 Wind stress
                !fkks = 0
                ks=2
                !fkif (lscdrag) ks = 1
                if (lscdrag) ks = 2

                do j=1,ngp
                    cdsdv     =  cds*denvvs(j,ks)
                    ustr(j,2) = -cdsdv*ua(j,kx)
                    vstr(j,2) = -cdsdv*va(j,kx)
                end do

            ! End of 'land-mode' computation
            end if

            ! Start of sea-sfc. heat fluxes computation
            ! 4.3 Sensible heat flux
            !fkks = 1
            ks=2
            chscp = chs*cp

            do j=1,ngp
                shf(j,2) = chscp*denvvs(j,ks)*(tsea(j)-t1(j,2))
            end do

            ! 4.4 Evaporation
            qsat0(:, 2) = q_sat(ngp, tsea, psa, rpe_literal(1.0))

            do j=1,ngp
                evap(j,2) = chs*denvvs(j,ks)*(qsat0(j,2)-q1(j,2))
            end do

            ! 4.5 Emission of lw radiation from the surface
            !     and net heat fluxes into sea surface
            do j=1,ngp
                slru(j,2)   = esbc*tsea(j)**4
                hfluxn(j,2) = ssrd(j)*(1.-alb_s(j))+slrd(j)-&
                    & (slru(j,2)+shf(j,2)+alhc*evap(j,2))
            end do

            ! End of sea-sfc. heat fluxes computation

            ! 3. Weighted average of surface fluxes and temperatures
            !    according to land-sea mask
            if (lfluxland)  then
                do j=1,ngp
                  ustr(j,3) = ustr(j,2)+fmask(j)*(ustr(j,1)-ustr(j,2))
                  vstr(j,3) = vstr(j,2)+fmask(j)*(vstr(j,1)-vstr(j,2))
                   shf(j,3) =  shf(j,2)+fmask(j)*( shf(j,1)- shf(j,2))
                  evap(j,3) = evap(j,2)+fmask(j)*(evap(j,1)-evap(j,2))
                  slru(j,3) = slru(j,2)+fmask(j)*(slru(j,1)-slru(j,2))
                end do

                do j=1,ngp
                  tsfc(j)  = tsea(j)+fmask(j)*(tland(j)-tsea(j))
                  tskin(j) = tsea(j)+fmask(j)*(tskin(j)-tsea(j))
                  t0(j)    = t1(j,2)+fmask(j)*(t1(j,1)- t1(j,2))
                  q0(j)    = q1(j,2)+fmask(j)*(q1(j,1)- q1(j,2))
                end do
            end if
        end

        subroutine sflset(phi0)
            ! subroutine sflset (phi0)
            !
            ! Purpose: compute orographic factor for land surface drag
            ! Input:   phi0   = surface geopotential            (2-dim)
            !          Initialized common blocks: sflfix

            use mod_physcon, only: gg
            use rp_emulator
            use mod_prec

            integer, parameter :: ngp=ix*il

            type(rpe_var), intent(in) :: phi0(ngp)
            integer :: j
            type(rpe_var) :: rhdrag

            rhdrag = 1./(gg*hdrag)

            do j=1,ngp
                forog(j)=1.+fhdrag*(1.-exp(-max(phi0(j),0.0_dp)*rhdrag))
            end do
        end subroutine sflset
end module surface_fluxes
