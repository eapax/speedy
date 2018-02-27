module mod_flx_sea
    use mod_atparam

    implicit none

    ! Fluxes at sea surface (all downward, except evaporation)
    ! Precipitation (sea)
    real, allocatable :: prec_s(:)
    
    ! Snowfall (sea)
    real, allocatable :: snowf_s(:)

    ! Evaporation (sea)
    real, allocatable :: evap_s(:)

    ! u-wind stress (sea)
    real, allocatable :: ustr_s(:)

    ! v-wind stress (sea)
    real, allocatable :: vstr_s(:)

    ! Sfc short-wave radiation (sea)
    real, allocatable :: ssr_s(:)

    ! Sfc long-wave radiation (sea)
    real, allocatable :: slr_s(:)

    ! Sensible heat flux (sea)
    real, allocatable :: shf_s(:)

    ! Latent heat flux (sea)
    real, allocatable :: ehf_s(:)

    ! Net heat flux into sea sfc.
    real, allocatable :: hflux_s(:)

    ! Net heat flux into sea-ice sfc.
    real, allocatable :: hflux_i(:)
    
    contains
        subroutine setup_flx_sea()
            allocate(prec_s(ix*il))    
            allocate(snowf_s(ix*il))        
            allocate(evap_s(ix*il))        
            allocate(ustr_s(ix*il))        
            allocate(vstr_s(ix*il))        
            allocate(ssr_s(ix*il))        
            allocate(slr_s(ix*il))        
            allocate(shf_s(ix*il))        
            allocate(ehf_s(ix*il))        
            allocate(hflux_s(ix*il))
            allocate(hflux_i(ix*il))
        end subroutine setup_flx_sea
end module
