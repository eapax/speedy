module mod_physvar
    use mod_atparam
    use rp_emulator

    implicit none

    ! ug1    = u-wind
    ! vg1    = v-wind
    ! tg1    = abs. temperature
    ! qg1    = specific humidity (g/kg)
    ! phig1  = geopotential
    ! pslg1  = log. of surface pressure
    type(rpe_var), dimension(:, :), allocatable :: ug1, vg1, tg1, qg1, phig1
    type(rpe_var), allocatable :: pslg1(:)

    ! se     = dry static energy
    ! rh     = relative humidity
    ! qsat   = saturation specific humidity (g/kg)
    type(rpe_var), dimension(:, :), allocatable :: se, rh, qsat

    ! psg    = surface pressure
    ! ts     = surface temperature
    ! tskin  = skin temperature
    ! u0     = near-surface u-wind
    ! v0     = near-surface v-wind
    ! t0     = near-surface air temperature
    ! q0     = near-surface specific humidity (g/kg)
    ! cloudc = total cloud cover (fraction)
    ! clstr  = stratiform cloud cover (fraction)
    ! cltop  = norm. pressure at cloud top
    ! prtop  = top of precipitation (level index)
    type(rpe_var), dimension(:), allocatable :: psg, ts, tskin, u0, v0, t0, q0, &
            cloudc, clstr, cltop, prtop

    ! tt_cnv  =  temperature tendency due to convection
    ! qt_cnv  = sp. humidity tendency due to convection
    ! tt_lsc  =  temperature tendency due to large-scale condensation
    ! qt_lsc  = sp. humidity tendency due to large-scale condensation
    ! tt_rsw  =  temperature tendency due to short-wave radiation
    ! tt_rlw  =  temperature tendency due to long-wave radiation
    ! ut_pbl  =       u-wind tendency due to PBL and diffusive processes
    ! vt_pbl  =       v-wind tendency due to PBL and diffusive processes
    ! tt_pbl  =  temperature tendency due to PBL and diffusive processes
    ! qt_pbl  = sp. humidity tendency due to PBL and diffusive processes
    type(rpe_var), dimension(:, :), allocatable :: tt_cnv, qt_cnv, tt_lsc, qt_lsc, &
            tt_rsw, tt_rlw, ut_pbl, vt_pbl, tt_pbl, qt_pbl

    ! precnv = convective precipitation  [g/(m^2 s)], total
    ! precls = large-scale precipitation [g/(m^2 s)], total
    ! snowcv = convective precipitation  [g/(m^2 s)], snow only
    ! snowls = large-scale precipitation [g/(m^2 s)], snow only
    ! cbmf   = cloud-base mass flux 
    ! tsr    = top-of-atm. shortwave radiation (downward)
    ! ssrd   = surface shortwave radiation (downward-only)
    ! ssr    = surface shortwave radiation (net downward)
    ! slrd   = surface longwave radiation  (downward-only)
    ! slr    = surface longwave radiation  (net upward) 
    ! olr    = outgoing longwave radiation (upward)
    ! slru   = surface longwave emission   (upward)
    !                                   (1:land, 2:sea, 3: wgt. average)
    ! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
    ! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
    ! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
    ! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
    ! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
    type(rpe_var), dimension(:), allocatable :: precnv, precls, snowcv, snowls, &
            cbmf, tsr, ssrd, ssr, slrd, slr, olr
    type(rpe_var), dimension(:, :), allocatable :: slru, ustr, vstr, shf, evap, hfluxn

    contains
        subroutine setup_physvar()
            allocate(ug1(ngp, kx))
            allocate(vg1(ngp, kx))
            allocate(tg1(ngp, kx))
            allocate(qg1(ngp, kx))
            allocate(phig1(ngp, kx))
            allocate(pslg1(ngp))
            allocate(se(ngp, kx))
            allocate(rh(ngp, kx))
            allocate(qsat(ngp, kx))
            allocate(psg(ngp))
            allocate(ts(ngp))
            allocate(tskin(ngp))
            allocate(u0(ngp))
            allocate(v0(ngp))
            allocate(t0(ngp))
            allocate(q0(ngp))
            allocate(cloudc(ngp))
            allocate(clstr(ngp))
            allocate(cltop(ngp))
            allocate(prtop(ngp))
            allocate(tt_cnv(ngp, kx))
            allocate(qt_cnv(ngp, kx))
            allocate(tt_lsc(ngp, kx))
            allocate(qt_lsc(ngp, kx))
            allocate(tt_rsw(ngp, kx))
            allocate(tt_rlw(ngp, kx))
            allocate(ut_pbl(ngp, kx))
            allocate(vt_pbl(ngp, kx))
            allocate(tt_pbl(ngp, kx))
            allocate(qt_pbl(ngp, kx))
            allocate(precnv(ngp))
            allocate(precls(ngp))
            allocate(snowcv(ngp))
            allocate(snowls(ngp))
            allocate(cbmf(ngp))
            allocate(tsr(ngp))
            allocate(ssrd(ngp))
            allocate(ssr(ngp))
            allocate(slrd(ngp))
            allocate(slr(ngp))
            allocate(olr(ngp))
            allocate(slru(ngp, 3))
            allocate(ustr(ngp, 3))
            allocate(vstr(ngp, 3))
            allocate(shf(ngp, 3))
            allocate(evap(ngp, 3))
            allocate(hfluxn(ngp, 3))
        end subroutine setup_physvar

        subroutine truncate_physvar()
            ug1 = ug1
            vg1 = vg1
            tg1 = tg1
            qg1 = qg1
            phig1 = phig1
            pslg1 = pslg1

            se = se
            rh = rh
            qsat = qsat

            psg = psg
            ts = ts
            tskin = tskin
            u0 = u0
            v0 = v0
            t0 = t0
            q0 = q0
            cloudc = cloudc
            clstr = clstr
            cltop = cltop
            prtop = prtop

            tt_cnv = tt_cnv
            qt_cnv = qt_cnv
            tt_lsc = tt_lsc
            qt_lsc = qt_lsc
            tt_rsw = tt_rsw
            tt_rlw = tt_rlw
            ut_pbl = ut_pbl
            vt_pbl = vt_pbl
            tt_pbl = tt_pbl
            qt_pbl = qt_pbl

            precnv = precnv
            precls = precls
            snowcv = snowcv
            snowls = snowls
            cbmf = cbmf
            tsr = tsr
            ssrd = ssrd
            ssr = ssr
            slrd = slrd
            slr = slr
            olr = olr
            slru = slru
            ustr = ustr
            vstr = vstr
            shf = shf
            evap = evap
            hfluxn = hfluxn
        end subroutine truncate_physvar
end module
