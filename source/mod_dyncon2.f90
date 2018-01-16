module mod_dyncon2
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public truncate_dyncon2
    public :: tref, tref1, tref2, tref3
    public :: xa, xb, xc, xd, xe, xf, xg, xh, xj, dhsx, elz

    ! Temp. profile for semi-imp. scheme (initial. in IMPINT)
    type(rpe_var), dimension(kx) :: tref, tref1, tref2, tref3

    type(rpe_var), dimension(kx,kx) :: xa, xb, xc, xd, xe
    type(rpe_var), dimension(kx,kx,lmax) :: xf, xg, xh, xj
    type(rpe_var) :: dhsx(kx), elz(mx,nx)

    contains
        subroutine truncate_dyncon2()
            tref = tref
            tref1 = tref1
            tref2 = tref2
            tref3 = tref3
            xa = xa
            xb = xb
            xc = xc
            xd = xd
            xe = xe
            xf = xf
            xg = xg
            xh = xh
            xj = xj
            dhsx = dhsx
            elz = elz
        end subroutine
end module
