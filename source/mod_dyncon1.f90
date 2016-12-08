module mod_dyncon1
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    private
    public rearth, omega, grav, akap, rgas, pi, a, g
    public hsg, dhs, fsg, dhsr, fsgr
    public radang, gsin, gcos, coriol
    public xgeop1, xgeop2

    ! Physical constants for dynamics
    real(dp), parameter :: rearth = 6.371e+6
    real(dp), parameter :: omega  = 7.292e-05
    real(dp), parameter :: grav   = 9.81
    real(dp), parameter :: akap   = 2./7.
    real(dp), parameter :: rgas   = akap*1004.
    real(dp), parameter :: pi = 4.*atan(1.)
    real(dp), parameter :: a  = rearth
    real(dp), parameter :: g  = grav

    ! Vertical level parameters (initial. in indyns)
    type(rpe_var) :: hsg(kxp), dhs(kx), fsg(kx), dhsr(kx), fsgr(kx)

    ! Functions of lat. and lon. (initial. in indyns)
    type(rpe_var) :: radang(il), gsin(il), gcos(il), coriol(il)

    ! Constants for hydrostatic eq. (initial. in indyns)
    type(rpe_var) :: xgeop1(kx), xgeop2(kx)
end module
