module mod_var_land
    use mod_atparam
    use rp_emulator

    implicit none

    ! Daily observed climatological fields over land
    type(rpe_var), allocatable :: stlcl_ob(:)     ! clim. land sfc. temperature
    type(rpe_var), allocatable :: snowdcl_ob(:)   ! clim. snow depth (water equiv)
    type(rpe_var), allocatable :: soilwcl_ob(:)   ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    type(rpe_var), allocatable :: stl_am(:)       ! land sfc. temperature
    type(rpe_var), allocatable :: snowd_am(:)     ! snow depth (water equiv)
    type(rpe_var), allocatable :: soilw_am(:)     ! soil water availability

    ! Land sfc. fields from land model
    type(rpe_var), allocatable :: stl_lm(:)       ! land-model sfc. temperature

    contains
        subroutine setup_var_land()
            allocate(stlcl_ob(ngp))
            allocate(snowdcl_ob(ngp))
            allocate(soilwcl_ob(ngp))
            allocate(stl_am(ngp))
            allocate(snowd_am(ngp))
            allocate(soilw_am(ngp))
            allocate(stl_lm(ngp))
        end subroutine setup_var_land

        subroutine truncate_var_land()
            stlcl_ob = stlcl_ob
            snowdcl_ob = snowdcl_ob
            soilwcl_ob = soilwcl_ob
            stl_am = stl_am
            snowd_am = snowd_am
            soilw_am = soilw_am
            stl_lm = stl_am
        end subroutine truncate_var_land
end module
