! Apply truncation to all constants used in speedy
subroutine truncate_rp()
    use mod_prec, only: set_precision

    use phy_convmf, only: truncate_convmf
    use phy_lscond, only: truncate_lscond
    use phy_cloud, only: truncate_cloud
    use phy_radsw, only: truncate_radsw
    use phy_radlw, only: truncate_radlw
    use phy_suflux, only: truncate_suflux
    use phy_vdifsc, only: truncate_vdifsc

    ! Truncate constants in individual physics schemes
    call set_precision('Convection')
    call truncate_convmf()
    call set_precision('Condensation')
    call truncate_lscond()
    call set_precision('Cloud')
    call truncate_cloud()
    call set_precision('Short-Wave Radiation')
    call truncate_radsw()
    call set_precision('Long-Wave Radiation')
    call truncate_radlw()
    call set_precision('Surface Fluxes')
    call truncate_suflux()
    call set_precision('Vertical Diffusion')
    call truncate_vdifsc()

    ! Set default precision to start model run
    call set_precision('Double')
end subroutine truncate_rp
