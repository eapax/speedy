! Copy douple precision coded parameters to rpe vars
subroutine ini_rp()
    use mod_dyncon1, only: init_dyncon1
    use mod_physcon, only: init_physcon

    call init_dyncon1
    call init_physcon
end subroutine ini_rp

! Apply truncation to all constants used in speedy
subroutine truncate_rp()
    use mod_cli_land, only: truncate_cli_land
    use mod_cli_sea, only: truncate_cli_sea
    use mod_cpl_land_model, only: truncate_land_model
    use mod_cplcon_sea, only: truncate_cplcon_sea
    use mod_cplvar_sea, only: truncate_cplvar_sea
    use mod_dyncon0, only: truncate_dyncon0
    use mod_dyncon1, only: truncate_dyncon1
    use mod_dyncon2, only: truncate_dyncon2
    use mod_fft, only: truncate_fft
    use mod_flx_land, only: truncate_flx_land
    use mod_flx_sea, only: truncate_flx_sea
    use mod_hdifcon, only: truncate_hdifcon
    use mod_physcon, only: truncate_physcon
    use mod_physvar, only: truncate_physvar
    use mod_randfor, only: truncate_randfor
    use mod_sppt, only: truncate_sppt
    use mod_surfcon, only: truncate_surfcon
    use mod_tsteps, only: truncate_tsteps
    use mod_var_land, only: truncate_var_land
    use mod_var_sea, only: truncate_var_sea
    use phy_convmf, only: truncate_convmf
    use phy_lscond, only: truncate_lscond
    use phy_radiat, only: truncate_radiat
    use phy_suflux, only: truncate_suflux
    use phy_vdifsc, only: truncate_vdifsc
    use spectral, only: truncate_spectral

    call truncate_cli_land()
    call truncate_cli_sea()
    call truncate_land_model()
    call truncate_cplcon_sea()
    call truncate_cplvar_sea()
    call truncate_dyncon0()
    call truncate_dyncon1()
    call truncate_dyncon2()
    call truncate_fft()
    call truncate_flx_land()
    call truncate_flx_sea()
    call truncate_hdifcon()
    call truncate_physcon()
    call truncate_physvar()
    call truncate_randfor()
    call truncate_sppt()
    call truncate_surfcon()
    call truncate_tsteps()
    call truncate_var_land()
    call truncate_var_sea()
    call truncate_convmf()
    call truncate_lscond()
    call truncate_radiat()
    call truncate_suflux()
    call truncate_vdifsc()
    call truncate_spectral()
end subroutine truncate_rp