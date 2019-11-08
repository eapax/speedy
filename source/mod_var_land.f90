module mod_var_land
    use mod_atparam
    use rp_emulator

    implicit none

    ! Daily observed climatological fields over land
    ! Interpolated from climatological fields in cpl_land.atm2land
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
            call apply_truncation(stlcl_ob)
            
            ! Denormal numbers at half precision but probably OK as zeros
            call apply_truncation(snowdcl_ob)

            call apply_truncation(soilwcl_ob)
            call apply_truncation(stl_am)
            
            ! Denormal numbers at half precision but probably OK as zeros
            call apply_truncation(snowd_am)
            
            call apply_truncation(soilw_am)
            call apply_truncation(stl_lm)
        end subroutine truncate_var_land
end module
