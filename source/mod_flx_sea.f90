module mod_flx_sea
    use mod_atparam

    implicit none

    ! Net heat flux into sea sfc.
    real, allocatable :: hflux_s(:)

    ! Net heat flux into sea-ice sfc.
    real, allocatable :: hflux_i(:)
    
    contains
        subroutine setup_flx_sea()
            allocate(hflux_s(ngp))
            allocate(hflux_i(ngp))
        end subroutine setup_flx_sea
end module
