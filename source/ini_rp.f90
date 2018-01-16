subroutine ini_rp
    use convection, only: init_cnvcon
    use condensation, only: init_lsccon
    use surface_fluxes, only: init_sflcon
    use vertical_diffusion, only: init_vdicon
    use mod_dyncon0, only: init_dyncon0
    use mod_dyncon1, only: init_dyncon1
    use mod_physcon, only: init_physcon
    use mod_radcon, only: init_radcon
    use mod_surfcon, only: init_surfcon
    use mod_tsteps, only: init_tsteps
    
    call init_cnvcon
    call init_dyncon0
    call init_dyncon1
    call init_lsccon
    call init_physcon
    call init_radcon
    call init_sflcon
    call init_surfcon
    call init_tsteps
    call init_vdicon
end subroutine

subroutine truncate_rp
    use mod_cli_land, only: truncate_cli_land
    use mod_cli_sea, only: truncate_cli_sea
    use mod_cpl_land_model, only: truncate_land_model
    use mod_cplcon_sea, only: truncate_cplcon_sea
    use mod_cplvar_sea, only: truncate_cplvar_sea
    use mod_dyncon1, only: truncate_dyncon1
    use mod_dyncon2, only: truncate_dyncon2
    use mod_fft, only: truncate_fft
    use mod_flx_land, only: truncate_flx_land
    use mod_flx_sea, only: truncate_flx_sea
    use mod_hdifcon, only: truncate_hdifcon
    use mod_physcon, only: truncate_physcon
    use mod_physvar, only: truncate_physvar
    use mod_radcon, only: truncate_radcon
    use mod_randfor, only: truncate_randfor
    use mod_spectral, only: truncate_spectral
    use mod_surfcon, only: truncate_surfcon
    use mod_var_land, only: truncate_var_land
    use mod_var_sea, only: truncate_var_sea
    use surface_fluxes, only: truncate_sflcon

    call truncate_cli_land
    call truncate_cli_sea
    call truncate_land_model
    call truncate_cplcon_sea
    call truncate_cplvar_sea
    call truncate_dyncon1
    call truncate_dyncon2
    call truncate_fft
    call truncate_flx_land
    call truncate_flx_sea
    call truncate_hdifcon
    call truncate_physcon
    call truncate_physvar
    call truncate_radcon
    call truncate_randfor
    call truncate_spectral
    call truncate_surfcon
    call truncate_var_land
    call truncate_var_sea
    call truncate_sflcon
end subroutine truncate_rp
