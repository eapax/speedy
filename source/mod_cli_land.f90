module mod_cli_land
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Land masks
    ! Fraction of land
    real(dp), allocatable :: fmask_l(:,:)

    ! Binary land mask
    real(dp), allocatable :: bmask_l(:,:)

    ! Monthly-mean climatological fields over land
    ! Loaded in ini_inbcon
    ! Used to interpolate to current date in cpl_land.atm2land
    ! Land surface temperature
    real(dp), allocatable :: stl12(:,:,:)

    ! Snow depth (water equiv.)
    real(dp), allocatable :: snowd12(:,:,:)

    ! Soil water availabilityend module
    real(dp), allocatable :: soilw12(:,:,:)

    contains
        subroutine setup_cli_land()
            allocate(fmask_l(ix, il))
            allocate(bmask_l(ix, il))
            allocate(stl12(ix, il, 12))
            allocate(snowd12(ix, il, 12))
            allocate(soilw12(ix, il, 12))
        end subroutine setup_cli_land
end module
