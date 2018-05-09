module mod_dyncon1
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    ! Physical constants for dynamics
    real(dp), parameter :: rearth_ = 6.371d+6
    real(dp), parameter :: omega_  = 7.292d-05
    real(dp), parameter :: grav_   = 9.81_dp
    real(dp), parameter :: akap_   = 2.0_dp/7.0_dp
    real(dp), parameter :: rgas_   = akap_*1004.0_dp
    real(dp), parameter :: pi_ = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: a_  = rearth_
    real(dp), parameter :: g_  = grav_

    ! Reduced precision versions
    type(rpe_var) :: rearth
    type(rpe_var) :: omega
    type(rpe_var) :: grav
    type(rpe_var) :: akap
    type(rpe_var) :: rgas
    type(rpe_var) :: pi
    type(rpe_var) :: a
    type(rpe_var) :: g

    ! Vertical level parameters (initial. in indyns)
    real, allocatable :: hsg(:), dhs(:), fsg(:), dhsr(:), fsgr(:)

    ! Functions of lat. and lon. (initial. in indyns)
    real, allocatable :: radang(:), gsin(:), gcos(:), coriol(:)

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
            allocate(gsin(il))
            allocate(gcos(il))
            allocate(coriol(il))
            allocate(xgeop1(kx))
            allocate(xgeop2(kx))
        end subroutine setup_dyncon1

        subroutine init_dyncon1()
            rearth = rearth_
            omega = omega_
            grav = grav_
            akap = akap_
            rgas = rgas_
            pi = pi_
            a = a_
            g = g_
        end subroutine init_dyncon1

        subroutine truncate_dyncon1()
            hsg = hsg
            dhs = dhs
            fsg = fsg
            dhsr = dhsr
            fsgr = fsgr
            radang = radang
            gsin = gsin
            gcos = gcos
            coriol = coriol
            xgeop1 = xgeop1
            xgeop2 = xgeop2
        end subroutine truncate_dyncon1

end module
