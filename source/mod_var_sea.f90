module mod_var_sea
    use mod_atparam
    use rp_emulator

    implicit none

    ! Daily observed climatological fields over sea
    ! Observed clim. SST
    type(rpe_var), allocatable :: sstcl_ob(:)

    ! Clim. sea ice fraction
    type(rpe_var), allocatable :: sicecl_ob(:)

    ! Clim. sea ice temperature
    type(rpe_var), allocatable :: ticecl_ob(:)

    ! Daily observed SST anomaly
    ! Observed SST anomaly
    type(rpe_var), allocatable :: sstan_ob(:)

    ! Daily climatological fields from ocean model
    ! Ocean model clim. SST
    type(rpe_var), allocatable :: sstcl_om(:)

    ! Sea sfc. fields used by atmospheric model
    ! SST (full-field)
    type(rpe_var), allocatable :: sst_am(:)

    ! SST anomaly
    type(rpe_var), allocatable :: sstan_am(:)

    ! Sea ice fraction
    type(rpe_var), allocatable :: sice_am(:)

    ! Sea ice temperature
    type(rpe_var), allocatable :: tice_am(:)

    ! Sea sfc. fields from ocean/sea-ice model
    ! Ocean model SST
    type(rpe_var), allocatable :: sst_om(:)

    ! Model sea ice fraction
    type(rpe_var), allocatable :: sice_om(:)

    ! Model sea ice temperature
    type(rpe_var), allocatable :: tice_om(:)

    ! Model SST + sea ice temp.
    type(rpe_var), allocatable :: ssti_om(:)

    ! Weight for obs. SST anomaly in coupled runs
    ! Weight mask for obs. SST
    type(rpe_var), allocatable :: wsst_ob(:)

    contains
        subroutine setup_var_sea()
            allocate(sstcl_ob(ngp))
            allocate(sicecl_ob(ngp))
            allocate(ticecl_ob(ngp))
            allocate(sstan_ob(ngp))
            allocate(sstcl_om(ngp))
            allocate(sst_am(ngp))
            allocate(sstan_am(ngp))
            allocate(sice_am(ngp))
            allocate(tice_am(ngp))
            allocate(sst_om(ngp))
            allocate(sice_om(ngp))
            allocate(tice_om(ngp))
            allocate(ssti_om(ngp))
            allocate(wsst_ob(ngp))
        end subroutine setup_var_sea

        subroutine truncate_var_sea()
            call apply_truncation(sstcl_ob)
            call apply_truncation(sicecl_ob)
            call apply_truncation(ticecl_ob)
            call apply_truncation(sstan_ob)
            call apply_truncation(sstcl_om)
            call apply_truncation(sst_am)
            call apply_truncation(sstan_am)
            call apply_truncation(sice_am)
            call apply_truncation(tice_am)
            call apply_truncation(sst_om)
            call apply_truncation(sice_om)
            call apply_truncation(tice_om)
            call apply_truncation(ssti_om)
            call apply_truncation(wsst_ob)
        end subroutine truncate_var_sea
end module
