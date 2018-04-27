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
            allocate(prec_s(ngp))    
            allocate(snowf_s(ngp))        
            allocate(evap_s(ngp))        
            allocate(ustr_s(ngp))        
            allocate(vstr_s(ngp))        
            allocate(ssr_s(ngp))        
            allocate(slr_s(ngp))        
            allocate(shf_s(ngp))        
            allocate(ehf_s(ngp))        
            allocate(hflux_s(ngp))
            allocate(hflux_i(ngp))
        end subroutine setup_flx_sea
end module
