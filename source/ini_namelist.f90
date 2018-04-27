subroutine ini_namelist()
    use mod_atparam, only: setup_resolution
    use downscaling, only: setup_downscaling
    use mod_cpl_flags, only: setup_coupling_flags
    use mod_cpl_land_model, only: setup_land
    use mod_cplcon_sea, only: setup_sea
    use mod_dyncon0, only: setup_dynamics
    use mod_physics, only: setup_physics
    use mod_surfcon, only: setup_surface
    use mod_tsteps, only: setup_timestepping
    use phy_convmf, only: setup_convection
    use phy_lscond, only: setup_condensation
    use phy_suflux, only: setup_surface_fluxes
    use phy_vdifsc, only: setup_vertical_diffusion
    use mod_sppt, only: setup_sppt

    use mod_cli_land, only: setup_cli_land
    use mod_cli_sea, only: setup_cli_sea
    use mod_cplvar_sea, only: setup_cplvar_sea
    use mod_dyncon1, only: setup_dyncon1
    use mod_dyncon2, only: setup_dyncon2
    use mod_dynvar, only: setup_dynvar
    use mod_fft, only: setup_fft
    use mod_flx_land, only: setup_flx_land
    use mod_flx_sea, only: setup_flx_sea
    use mod_hdifcon, only: setup_hdifcon
    use mod_physcon, only: setup_physcon
    use mod_physvar, only: setup_physvar
    use mod_radcon, only: setup_radcon
    use mod_randfor, only: setup_randfor
    use mod_var_land, only: setup_var_land
    use mod_var_sea, only: setup_var_sea
    use spectral, only: setup_spectral

    open(99, file='speedy.nml')

    ! Read parameters from the namelist and allocate arrays
    ! Call setup_resolution first as it contains the dimensions needed for
    ! allocating arrays
    call setup_resolution(99)
    call setup_downscaling(99)
    call setup_coupling_flags(99)
    call setup_land(99)
    call setup_sea(99)
    call setup_dynamics(99)
    call setup_physics(99)
    call setup_surface(99)
    call setup_timestepping(99)
    call setup_convection(99)
    call setup_condensation(99)
    call setup_surface_fluxes(99)
    call setup_vertical_diffusion(99)
    call setup_sppt(99)

    close(99)

    print *,'Read namelist'

    ! Remaining setup routines only allocate arrays
    call setup_cli_land()
    call setup_cli_sea()
    call setup_cplvar_sea()
    call setup_dyncon1()
    call setup_dyncon2()
    call setup_dynvar()
    call setup_fft()
    call setup_flx_land()
    call setup_flx_sea()
    call setup_hdifcon()
    call setup_physcon()
    call setup_physvar()
    call setup_radcon()
    call setup_randfor()
    call setup_var_land()
    call setup_var_sea()
    call setup_spectral()
end subroutine ini_namelist