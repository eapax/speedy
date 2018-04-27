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
            allocate(prec_l(ngp))
            allocate(snowf_l(ngp))
            allocate(evap_l(ngp))
            allocate(ustr_l(ngp))
            allocate(vstr_l(ngp))
            allocate(ssr_l(ngp))
            allocate(slr_l(ngp))
            allocate(shf_l(ngp))
            allocate(ehf_l(ngp))
            allocate(hflux_l(ngp))
        end subroutine setup_flx_land
end module
