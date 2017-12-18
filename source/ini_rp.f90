subroutine ini_rp
    use convection, only: init_cnvcon
    use condensation, only: init_lsccon
    use surface_fluxes, only: init_sflcon
    use vertical_diffusion, only: init_vdicon
    use mod_dyncon0
    use mod_dyncon1
    use mod_physcon
    use mod_radcon
    use mod_surfcon
    use mod_tsteps
    
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
