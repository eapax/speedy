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
            sstcl_ob = sstcl_ob
            sicecl_ob = sicecl_ob
            ticecl_ob = ticecl_ob
            sstan_ob = sstan_ob
            sstcl_om = sstcl_om
            sst_am = sst_am
            sstan_am = sstan_am
            sice_am = sice_am
            tice_am = tice_am
            sst_om = sst_om
            sice_om = sice_om
            tice_om = tice_om
            ssti_om = ssti_om
            wsst_ob = wsst_ob
        end subroutine truncate_var_sea
end module
