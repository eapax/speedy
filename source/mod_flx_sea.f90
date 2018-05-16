module mod_flx_sea
    use mod_atparam
	use rp_emulator

    implicit none

    ! Net heat flux into sea sfc.
    type(rpe_var), allocatable :: hflux_s(:)

    ! Net heat flux into sea-ice sfc.
    type(rpe_var), allocatable :: hflux_i(:)

    contains
        subroutine setup_flx_sea()
            allocate(hflux_s(ngp))
            allocate(hflux_i(ngp))
        end subroutine setup_flx_sea

        subroutine truncate_flx_sea()
            call apply_truncation(hflux_s)
            call apply_truncation(hflux_i)
        end subroutine truncate_flx_sea
end module
