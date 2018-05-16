module mod_flx_land
    use mod_atparam
    use rp_emulator

    implicit none

    ! Fluxes at land surface (all downward, except evaporation)
    ! Precipitation (land)
    type(rpe_var), allocatable :: prec_l(:)

    ! Snowfall (land)
    type(rpe_var), allocatable :: snowf_l(:)

    ! Evaporation (land)
    type(rpe_var), allocatable :: evap_l(:)

    ! u-wind stress (land)
    type(rpe_var), allocatable :: ustr_l(:)

    ! v-wind stress (land)
    type(rpe_var), allocatable :: vstr_l(:)

    ! Sfc short-wave radiation (land)
    type(rpe_var), allocatable :: ssr_l(:)

    ! Sfc long-wave radiation (land)
    type(rpe_var), allocatable :: slr_l(:)

    ! Sensible heat flux (land)
    type(rpe_var), allocatable :: shf_l(:)

    ! Latent heat flux (land)
    type(rpe_var), allocatable :: ehf_l(:)

    ! Net heat flux into land sfc.end module
    type(rpe_var), allocatable :: hflux_l(:)

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

        subroutine truncate_flx_land()
            call apply_truncation(prec_l)
            call apply_truncation(snowf_l)
            call apply_truncation(evap_l)
            call apply_truncation(ustr_l)
            call apply_truncation(vstr_l)
            call apply_truncation(ssr_l)
            call apply_truncation(slr_l)
            call apply_truncation(shf_l)
            call apply_truncation(ehf_l)
            call apply_truncation(hflux_l)
        end subroutine truncate_flx_land
end module
