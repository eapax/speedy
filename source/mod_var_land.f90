module mod_var_land
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public stlcl_ob, snowdcl_ob, soilwcl_ob, stl_am, snowd_am, soilw_am, stl_lm
    public truncate_var_land

    ! Daily observed climatological fields over land
    type(rpe_var) :: stlcl_ob(ix*il)              ! clim. land sfc. temperature 
    type(rpe_var) :: snowdcl_ob(ix*il)              ! clim. snow depth (water equiv)
    type(rpe_var) :: soilwcl_ob(ix*il)              ! clim. soil water availability

    ! Land sfc. fields used by atmospheric model
    type(rpe_var) :: stl_am(ix*il)                 ! land sfc. temperature
    type(rpe_var) :: snowd_am(ix*il)                 ! snow depth (water equiv)
    type(rpe_var) :: soilw_am(ix*il)                 ! soil water availability

    ! Land sfc. fields from land model
    type(rpe_var) :: stl_lm(ix*il)                 ! land-model sfc. temperature

    contains
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
