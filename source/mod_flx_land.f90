module mod_flx_land
    use mod_atparam

    implicit none

    ! Net heat flux into land sfc.end module
    real, allocatable :: hflux_l(:)
    
    contains
        subroutine setup_flx_land()
            allocate(hflux_l(ngp))
        end subroutine setup_flx_land
end module
