module mod_flx_land
    use mod_atparam
    use rp_emulator

    implicit none

    ! Net heat flux into land sfc.end module
    type(rpe_var), allocatable :: hflux_l(:)

    contains
        subroutine setup_flx_land()
            allocate(hflux_l(ngp))
        end subroutine setup_flx_land

        subroutine truncate_flx_land()
            call apply_truncation(hflux_l)
        end subroutine truncate_flx_land
end module
