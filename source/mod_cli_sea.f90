module mod_cli_sea
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Sea masks
    ! Fraction of sea
    real(dp), allocatable :: fmask_s(:,:)

    ! Binary sea mask
    real(dp), allocatable :: bmask_s(:,:)

    ! Grid latitudes
    real(dp), allocatable :: deglat_s(:)

    ! Monthly-mean climatological fields over sea
    ! Sea/ice surface temperature
    real(dp), allocatable :: sst12(:,:,:)

    ! Sea ice fraction
    real(dp), allocatable :: sice12(:,:,:)

    ! SST anomaly fields
    ! SST anomaly in 3 consecutive months
    real(dp), allocatable :: sstan3(:,:,:)

    ! Climatological fields from model output
    ! Annual-mean heat flux into sea sfc.
    real(dp), allocatable :: hfseacl(:,:)

    ! Ocean model SST climatology
    real(dp), allocatable :: sstom12(:,:,:)

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
