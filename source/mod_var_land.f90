module mod_var_land
    use mod_atparam

    implicit none

    ! Daily observed climatological fields over land
    real, allocatable :: stlcl_ob(:)     ! clim. land sfc. temperature
    real, allocatable :: snowdcl_ob(:)   ! clim. snow depth (water equiv)
    real, allocatable :: soilwcl_ob(:)   ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real, allocatable :: stl_am(:)       ! land sfc. temperature
    real, allocatable :: snowd_am(:)     ! snow depth (water equiv)
    real, allocatable :: soilw_am(:)     ! soil water availability

    ! Land sfc. fields from land model
    real, allocatable :: stl_lm(:)       ! land-model sfc. temperature
    
    contains
        subroutine setup_var_land()
            allocate(stlcl_ob(ix*il))
            allocate(snowdcl_ob(ix*il))
            allocate(soilwcl_ob(ix*il))
            allocate(stl_am(ix*il))
            allocate(snowd_am(ix*il))
            allocate(soilw_am(ix*il))
            allocate(stl_lm(ix*il))
        end subroutine setup_var_land
end module
