module mod_var_sea
    use mod_atparam

    implicit none

    ! Daily observed climatological fields over sea
    ! Observed clim. SST
    real, allocatable :: sstcl_ob(:)

    ! Clim. sea ice fraction
    real, allocatable :: sicecl_ob(:)

    ! Clim. sea ice temperature
    real, allocatable :: ticecl_ob(:)

    ! Daily observed SST anomaly
    ! Observed SST anomaly
    real, allocatable :: sstan_ob(:)

    ! Daily climatological fields from ocean model
    ! Ocean model clim. SST
    real, allocatable :: sstcl_om(:)

    ! Sea sfc. fields used by atmospheric model
    ! SST (full-field)
    real, allocatable :: sst_am(:)

    ! SST anomaly
    real, allocatable :: sstan_am(:)

    ! Sea ice fraction
    real, allocatable :: sice_am(:)

    ! Sea ice temperature
    real, allocatable :: tice_am(:)

    ! Sea sfc. fields from ocean/sea-ice model
    ! Ocean model SST
    real, allocatable :: sst_om(:)

    ! Model sea ice fraction
    real, allocatable :: sice_om(:)

    ! Model sea ice temperature
    real, allocatable :: tice_om(:)

    ! Model SST + sea ice temp.
    real, allocatable :: ssti_om(:)

    ! Weight for obs. SST anomaly in coupled runs
    ! Weight mask for obs. SST
    real, allocatable :: wsst_ob(:)
    
    contains
        subroutine setup_var_sea()
            allocate(sstcl_ob(ix*il))
            allocate(sicecl_ob(ix*il))
            allocate(ticecl_ob(ix*il))
            allocate(sstan_ob(ix*il))
            allocate(sstcl_om(ix*il))
            allocate(sst_am(ix*il))
            allocate(sstan_am(ix*il))
            allocate(sice_am(ix*il))
            allocate(tice_am(ix*il))
            allocate(sst_om(ix*il))
            allocate(sice_om(ix*il))
            allocate(tice_om(ix*il))
            allocate(ssti_om(ix*il))
            allocate(wsst_ob(ix*il))
        end subroutine setup_var_sea
end module
