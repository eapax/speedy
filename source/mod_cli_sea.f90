module mod_cli_sea
    use mod_atparam

    implicit none

    ! Sea masks
    ! Fraction of sea
    real, allocatable :: fmask_s(:,:)

    ! Binary sea mask
    real, allocatable :: bmask_s(:,:)

    ! Grid latitudes
    real, allocatable :: deglat_s(:)

    ! Monthly-mean climatological fields over sea
    ! Sea/ice surface temperature
    real, allocatable :: sst12(:,:,:)

    ! Sea ice fraction
    real, allocatable :: sice12(:,:,:)

    ! SST anomaly fields
    ! SST anomaly in 3 consecutive months
    real, allocatable :: sstan3(:,:,:)

    ! Climatological fields from model output
    ! Annual-mean heat flux into sea sfc.
    real, allocatable :: hfseacl(:,:)

    ! Ocean model SST climatology
    real, allocatable :: sstom12(:,:,:)
    
    contains
        subroutine setup_cli_sea()
            allocate(fmask_s(ix, il))
            allocate(bmask_s(ix, il))
            allocate(deglat_s(il))
            allocate(sst12(ix, il, 12))
            allocate(sice12(ix, il, 12))
            allocate(sstan3(ix, il, 3))
            allocate(hfseacl(ix, il))
            allocate(sstom12(ix, il, 12))
        end subroutine setup_cli_sea
end module
