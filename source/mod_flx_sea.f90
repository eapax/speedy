module mod_flx_sea
    use mod_atparam
	use rp_emulator

    implicit none

    ! Fluxes at sea surface (all downward, except evaporation)
    ! Precipitation (sea)
    type(rpe_var), allocatable :: prec_s(:)
    
    ! Snowfall (sea)
    type(rpe_var), allocatable :: snowf_s(:)

    ! Evaporation (sea)
    type(rpe_var), allocatable :: evap_s(:)

    ! u-wind stress (sea)
    type(rpe_var), allocatable :: ustr_s(:)

    ! v-wind stress (sea)
    type(rpe_var), allocatable :: vstr_s(:)

    ! Sfc short-wave radiation (sea)
    type(rpe_var), allocatable :: ssr_s(:)

    ! Sfc long-wave radiation (sea)
    type(rpe_var), allocatable :: slr_s(:)

    ! Sensible heat flux (sea)
    type(rpe_var), allocatable :: shf_s(:)

    ! Latent heat flux (sea)
    type(rpe_var), allocatable :: ehf_s(:)

    ! Net heat flux into sea sfc.
    type(rpe_var), allocatable :: hflux_s(:)

    ! Net heat flux into sea-ice sfc.
    type(rpe_var), allocatable :: hflux_i(:)

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

        subroutine truncate_flx_sea()
            call apply_truncation(prec_s)
            call apply_truncation(snowf_s)
            call apply_truncation(evap_s)
            call apply_truncation(ustr_s)
            call apply_truncation(vstr_s)
            call apply_truncation(ssr_s)
            call apply_truncation(slr_s)
            call apply_truncation(shf_s)
            call apply_truncation(ehf_s)
            call apply_truncation(hflux_s)
            call apply_truncation(hflux_i)
        end subroutine truncate_flx_sea
end module
