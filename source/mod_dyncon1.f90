module mod_dyncon1
    use mod_atparam

    implicit none

    ! Physical constants for dynamics
    real, parameter :: rearth = 6.371e+6
    real, parameter :: omega  = 7.292e-05
    real, parameter :: grav   = 9.81
    real, parameter :: akap   = 2./7.
    real, parameter :: rgas   = akap*1004.
    real, parameter :: pi = 4.*atan(1.)
    real, parameter :: a  = rearth
    real, parameter :: g  = grav

    ! Vertical level parameters (initial. in indyns)
    real, allocatable :: hsg(:), dhs(:), fsg(:), dhsr(:), fsgr(:)

    ! Functions of lat. and lon. (initial. in indyns)
    real, allocatable :: radang(:), coriol(:)

    ! Constants for hydrostatic eq. (initial. in indyns)
    real, allocatable :: xgeop1(:), xgeop2(:)
    
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
