module mod_var_land
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Daily observed climatological fields over land
    ! Interpolated from climatological fields in cpl_land.atm2land
    real(dp), allocatable :: stlcl_ob(:)     ! clim. land sfc. temperature
    real(dp), allocatable :: snowdcl_ob(:)   ! clim. snow depth (water equiv)
    real(dp), allocatable :: soilwcl_ob(:)   ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    real(dp), allocatable :: stl_am(:)       ! land sfc. temperature
    real(dp), allocatable :: snowd_am(:)     ! snow depth (water equiv)
    real(dp), allocatable :: soilw_am(:)     ! soil water availability

    ! Land sfc. fields from land model
    real(dp), allocatable :: stl_lm(:)       ! land-model sfc. temperature

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
end module
