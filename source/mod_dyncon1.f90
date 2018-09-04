module mod_dyncon1
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Physical constants for dynamics setup
    real(dp), parameter :: rearth = 6.371d+6
    real(dp), parameter :: omega  = 7.292d-05
    real(dp), parameter :: grav   = 9.81_dp

    ! Physical constants for dynamics integration and setup
    real(dp), parameter :: akap   = 2.0_dp/7.0_dp
    real(dp), parameter :: rgas   = akap*1004.0_dp

    ! Vertical level parameters (initial. in indyns)
    real(dp), allocatable :: hsg(:), dhs(:), fsg(:), dhsr(:), fsgr(:)

    ! Functions of lat. and lon. (initial. in indyns)
    real(dp), allocatable :: radang(:), coriol(:)

    ! Constants for hydrostatic eq. (initial. in indyns)
    real(dp), allocatable :: xgeop1(:), xgeop2(:)

    contains
        subroutine setup_dyncon1()
            allocate(hsg(kxp))
            allocate(dhs(kx))
            allocate(fsg(kx))
            allocate(dhsr(kx))
            allocate(fsgr(kx))
            allocate(radang(il))
            allocate(coriol(il))
            allocate(xgeop1(kx))
            allocate(xgeop2(kx))
        end subroutine setup_dyncon1
end module
