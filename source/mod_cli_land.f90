module mod_cli_land
    use mod_atparam
    use rp_emulator

    implicit none

    ! Land masks
    ! Fraction of land
    type(rpe_var), allocatable :: fmask_l(:,:)

    ! Binary land mask
    type(rpe_var), allocatable :: bmask_l(:,:)

    ! Monthly-mean climatological fields over land
    ! Land surface temperature
    type(rpe_var), allocatable :: stl12(:,:,:)

    ! Snow depth (water equiv.)
    type(rpe_var), allocatable :: snowd12(:,:,:)

    ! Soil water availabilityend module
    type(rpe_var), allocatable :: soilw12(:,:,:)

    contains
        subroutine setup_cli_land()
            allocate(fmask_l(ix, il))
            allocate(bmask_l(ix, il))
            allocate(stl12(ix, il, 12))
            allocate(snowd12(ix, il, 12))
            allocate(soilw12(ix, il, 12))
        end subroutine setup_cli_land

        subroutine truncate_cli_land()
            fmask_l = fmask_l
            bmask_l = bmask_l
            stl12 = stl12
            snowd12 = snowd12
            soilw12 = soilw12
        end subroutine truncate_cli_land
end module
