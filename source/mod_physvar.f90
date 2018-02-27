module mod_physvar
    use mod_atparam

    implicit none

    ! ug1    = u-wind
    ! vg1    = v-wind
    ! tg1    = abs. temperature
    ! qg1    = specific humidity (g/kg)
    ! phig1  = geopotential
    ! pslg1  = log. of surface pressure
    real, dimension(:, :), allocatable :: ug1, vg1, tg1, qg1, phig1
    real, allocatable :: pslg1(:)

    ! se     = dry static energy
    ! rh     = relative humidity
    ! qsat   = saturation specific humidity (g/kg)
    real, dimension(:, :), allocatable :: se, rh, qsat

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
    real, dimension(:), allocatable :: psg, ts, tskin, u0, v0, t0, q0, &
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
    real, dimension(:, :), allocatable :: tt_cnv, qt_cnv, tt_lsc, qt_lsc, &
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
    real, dimension(:), allocatable :: precnv, precls, snowcv, snowls, &
            cbmf, tsr, ssrd, ssr, slrd, slr, olr
    real, dimension(:, :), allocatable :: slru, ustr, vstr, shf, evap, hfluxn
    
    contains
        subroutine setup_physvar()
            allocate(ug1(ix*il, kx))
            allocate(vg1(ix*il, kx))
            allocate(tg1(ix*il, kx))
            allocate(qg1(ix*il, kx))
            allocate(phig1(ix*il, kx))
            allocate(pslg1(ix*il))
            allocate(se(ix*il, kx))
            allocate(rh(ix*il, kx))
            allocate(qsat(ix*il, kx))
            allocate(psg(ix*il))
            allocate(ts(ix*il))
            allocate(tskin(ix*il))
            allocate(u0(ix*il))
            allocate(v0(ix*il))
            allocate(t0(ix*il))
            allocate(q0(ix*il))
            allocate(cloudc(ix*il))
            allocate(clstr(ix*il))
            allocate(cltop(ix*il))
            allocate(prtop(ix*il))
            allocate(tt_cnv(ix*il, kx))
            allocate(qt_cnv(ix*il, kx))
            allocate(tt_lsc(ix*il, kx))
            allocate(qt_lsc(ix*il, kx))
            allocate(tt_rsw(ix*il, kx))
            allocate(tt_rlw(ix*il, kx))
            allocate(ut_pbl(ix*il, kx))
            allocate(vt_pbl(ix*il, kx))
            allocate(tt_pbl(ix*il, kx))
            allocate(qt_pbl(ix*il, kx))
            allocate(precnv(ix*il))
            allocate(precls(ix*il))
            allocate(snowcv(ix*il))
            allocate(snowls(ix*il))
            allocate(cbmf(ix*il))
            allocate(tsr(ix*il))
            allocate(ssrd(ix*il))
            allocate(ssr(ix*il))
            allocate(slrd(ix*il))
            allocate(slr(ix*il))
            allocate(olr(ix*il))
            allocate(slru(ix*il, 3))
            allocate(ustr(ix*il, 3))
            allocate(vstr(ix*il, 3))
            allocate(shf(ix*il, 3))
            allocate(evap(ix*il, 3))
            allocate(hfluxn(ix*il, 3))
        end subroutine setup_physvar
end module
