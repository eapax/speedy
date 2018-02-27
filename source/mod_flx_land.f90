module mod_flx_land
    use mod_atparam

    implicit none

    ! Fluxes at land surface (all downward, except evaporation)
    ! Precipitation (land)
    real, allocatable :: prec_l(:)

    ! Snowfall (land)
    real, allocatable :: snowf_l(:)

    ! Evaporation (land)
    real, allocatable :: evap_l(:)

    ! u-wind stress (land)
    real, allocatable :: ustr_l(:)

    ! v-wind stress (land)
    real, allocatable :: vstr_l(:)

    ! Sfc short-wave radiation (land)
    real, allocatable :: ssr_l(:)

    ! Sfc long-wave radiation (land)
    real, allocatable :: slr_l(:)

    ! Sensible heat flux (land)
    real, allocatable :: shf_l(:)

    ! Latent heat flux (land)
    real, allocatable :: ehf_l(:)

    ! Net heat flux into land sfc.end module
    real, allocatable :: hflux_l(:)
    
    contains
        subroutine setup_flx_land()
            allocate(prec_l(ix*il))        
            allocate(snowf_l(ix*il))        
            allocate(evap_l(ix*il))        
            allocate(ustr_l(ix*il))        
            allocate(vstr_l(ix*il))        
            allocate(ssr_l(ix*il))        
            allocate(slr_l(ix*il))        
            allocate(shf_l(ix*il))        
            allocate(ehf_l(ix*il))        
            allocate(hflux_l(ix*il))        
        end subroutine setup_flx_land
end module
