module mod_cli_sea
    use mod_atparam
    use rp_emulator

    implicit none

    ! Sea masks
    ! Fraction of sea
    type(rpe_var), allocatable :: fmask_s(:,:)

    ! Binary sea mask
    type(rpe_var), allocatable :: bmask_s(:,:)

    ! Grid latitudes
    type(rpe_var), allocatable :: deglat_s(:)

    ! Monthly-mean climatological fields over sea
    ! Sea/ice surface temperature
    type(rpe_var), allocatable :: sst12(:,:,:)

    ! Sea ice fraction
    type(rpe_var), allocatable :: sice12(:,:,:)

    ! SST anomaly fields
    ! SST anomaly in 3 consecutive months
    type(rpe_var), allocatable :: sstan3(:,:,:)

    ! Climatological fields from model output
    ! Annual-mean heat flux into sea sfc.
    type(rpe_var), allocatable :: hfseacl(:,:)

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

        subroutine truncate_cli_sea()
            fmask_s = fmask_s
            bmask_s = bmask_s
            deglat_s = deglat_s
            sst12 = sst12
            sice12 = sice12
            sstan3 = sstan3
            hfseacl = hfseacl
            sstom12 = sstom12
        end subroutine truncate_cli_sea
end module
