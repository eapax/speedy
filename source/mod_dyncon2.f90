module mod_dyncon2
    use mod_atparam
    use rp_emulator

    implicit none

    ! Temp. profile for semi-imp. scheme (initial. in IMPINT)
    type(rpe_var), dimension(:),       allocatable :: tref, tref1, tref2, tref3
    type(rpe_var), dimension(:, :),    allocatable :: xa, xb, xc, xd, xe
    type(rpe_var), dimension(:, :, :), allocatable :: xf, xg, xh, xj
    type(rpe_var), allocatable :: dhsx(:), elz(:, :)

    contains
        subroutine setup_dyncon2()
            allocate(tref(kx))
            allocate(tref1(kx))
            allocate(tref2(kx))
            allocate(tref3(kx))

            allocate(xa(kx, kx))
            allocate(xb(kx, kx))
            allocate(xc(kx, kx))
            allocate(xd(kx, kx))
            allocate(xe(kx, kx))

            allocate(xf(kx, kx, lmax))
            allocate(xg(kx, kx, lmax))
            allocate(xh(kx, kx, lmax))
            allocate(xj(kx, kx, lmax))

            allocate(dhsx(kx))
            allocate(elz(mx,nx))
        end subroutine setup_dyncon2

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
        end subroutine truncate_dyncon2
end module
