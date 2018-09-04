module mod_var_sea
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Daily observed climatological fields over sea
    ! Observed clim. SST
    real(dp), allocatable :: sstcl_ob(:)

    ! Clim. sea ice fraction
    real(dp), allocatable :: sicecl_ob(:)

    ! Clim. sea ice temperature
    real(dp), allocatable :: ticecl_ob(:)

    ! Daily observed SST anomaly
    ! Observed SST anomaly
    real(dp), allocatable :: sstan_ob(:)

    ! Daily climatological fields from ocean model
    ! Ocean model clim. SST
    real(dp), allocatable :: sstcl_om(:)

    ! Sea sfc. fields used by atmospheric model
    ! SST (full-field)
    real(dp), allocatable :: sst_am(:)

    ! SST anomaly
    real(dp), allocatable :: sstan_am(:)

    ! Sea ice fraction
    real(dp), allocatable :: sice_am(:)

    ! Sea ice temperature
    real(dp), allocatable :: tice_am(:)

    ! Sea sfc. fields from ocean/sea-ice model
    ! Ocean model SST
    real(dp), allocatable :: sst_om(:)

    ! Model sea ice fraction
    real(dp), allocatable :: sice_om(:)

    ! Model sea ice temperature
    real(dp), allocatable :: tice_om(:)

    ! Model SST + sea ice temp.
    real(dp), allocatable :: ssti_om(:)

    ! Weight for obs. SST anomaly in coupled runs
    ! Weight mask for obs. SST
    real(dp), allocatable :: wsst_ob(:)

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
end module
