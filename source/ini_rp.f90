! Copy douple precision coded parameters to rpe vars
subroutine ini_rp()
    use mod_dyncon1, only: init_dyncon1
    use mod_physcon, only: init_physcon

    call init_dyncon1
    call init_physcon
end subroutine ini_rp

! Apply truncation to all constants used in speedy
subroutine truncate_rp()
    use mod_prec, only: set_precision
    use mod_cli_land, only: truncate_cli_land
    use mod_cli_sea, only: truncate_cli_sea
    use mod_cpl_land_model, only: truncate_land_model
    use mod_cplcon_sea, only: truncate_cplcon_sea
    use mod_cplvar_sea, only: truncate_cplvar_sea
    use mod_dyncon0, only: truncate_dyncon0
    use mod_dyncon1, only: truncate_dyncon1
    use mod_dyncon2, only: truncate_dyncon2
    use mod_fft, only: truncate_fft
    use mod_fordate, only: truncate_fordate
    use mod_hdifcon, only: truncate_hdifcon
    use mod_physcon, only: truncate_physcon
    use mod_solar, only: truncate_solar
    use mod_surfcon, only: truncate_surfcon
    use mod_tsteps, only: truncate_tsteps
    use mod_var_land, only: truncate_var_land
    use mod_var_sea, only: truncate_var_sea

    use phy_convmf, only: truncate_convmf
    use phy_lscond, only: truncate_lscond
    use phy_cloud, only: truncate_cloud
    use phy_radsw, only: truncate_radsw
    use phy_radlw, only: truncate_radlw
    use phy_suflux, only: truncate_suflux
    use phy_vdifsc, only: truncate_vdifsc
    use phy_sppt, only: truncate_sppt

    use spectral, only: truncate_spectral

    call set_precision('Default')
    ! Truncate climatological fields used in surface model
    call truncate_cli_land()
    call truncate_cli_sea()
    ! Truncate constants and variables in land model
    call truncate_land_model()
    call truncate_var_land()
    ! Truncate constants in sea model
    call truncate_cplcon_sea()
    ! Truncate variables in sea model
    call truncate_cplvar_sea()
    call truncate_var_sea()
    ! Truncate constants for determining forcing
    call truncate_fordate()
    call truncate_solar()

    ! Truncate general dynamics constants used in multiple schemes
    call set_precision('Default')
    call truncate_dyncon0()
    call truncate_dyncon1()
    call truncate_dyncon2()
    call truncate_hdifcon()

    ! Truncate timestepping constants
    call set_precision('Timestepping')
    call truncate_tsteps()

    ! Truncate FFT wsave array
    call set_precision('Spectral Transform')
    call truncate_fft()
    ! Truncate spectral transform constants
    call truncate_spectral()

    ! Truncate general physics constants used in multiple schemes
    call set_precision('Grid Physics')
    call truncate_physcon()
    call truncate_surfcon()

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
    call set_precision('SPPT')
    call truncate_sppt()

    ! Set default precision to start model run
    call set_precision('Default')
end subroutine truncate_rp
