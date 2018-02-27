module mod_cli_land
    use mod_atparam

    implicit none

    ! Land masks
    ! Fraction of land
    real, allocatable :: fmask_l(:,:)

    ! Binary land mask
    real, allocatable :: bmask_l(:,:)

    ! Monthly-mean climatological fields over land
    ! Land surface temperature
    real, allocatable :: stl12(:,:,:)

    ! Snow depth (water equiv.)
    real, allocatable :: snowd12(:,:,:)

    ! Soil water availabilityend module
    real, allocatable :: soilw12(:,:,:)

    contains
        subroutine setup_cli_land()
            allocate(fmask_l(ix, il))
            allocate(bmask_l(ix, il))
            allocate(stl12(ix, il, 12))
            allocate(snowd12(ix, il, 12))
            allocate(soilw12(ix, il, 12))
        end subroutine setup_cli_land
end module
