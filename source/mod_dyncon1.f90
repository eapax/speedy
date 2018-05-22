module mod_dyncon1
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    ! Physical constants for dynamics setup
    real(dp), parameter :: rearth = 6.371d+6
    real(dp), parameter :: omega  = 7.292d-05
    real(dp), parameter :: grav   = 9.81_dp

    ! Physical constants for dynamics integration and setup
    real(dp), parameter :: akap_   = 2.0_dp/7.0_dp
    real(dp), parameter :: rgas_   = akap_*1004.0_dp
    ! Reduced precision versions
    type(rpe_var) :: akap
    type(rpe_var) :: rgas

    ! Vertical level parameters (initial. in indyns)
    type(rpe_var), allocatable :: hsg(:), dhs(:), fsg(:), dhsr(:), fsgr(:)

    ! Functions of lat. and lon. (initial. in indyns)
    type(rpe_var), allocatable :: radang(:), coriol(:)

    ! Constants for hydrostatic eq. (initial. in indyns)
    type(rpe_var), allocatable :: xgeop1(:), xgeop2(:)

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

        subroutine init_dyncon1()
            akap = akap_
            rgas = rgas_
        end subroutine init_dyncon1

        subroutine truncate_dyncon1()
            call apply_truncation(akap)
            call apply_truncation(rgas)
            call apply_truncation(hsg)
            call apply_truncation(dhs)
            call apply_truncation(fsg)
            call apply_truncation(dhsr)
            call apply_truncation(fsgr)
            call apply_truncation(radang)
            call apply_truncation(coriol)
            call apply_truncation(xgeop1)
            call apply_truncation(xgeop2)
        end subroutine truncate_dyncon1

end module
