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
            call apply_truncation(tref)
            call apply_truncation(tref1)
            call apply_truncation(tref2)
            call apply_truncation(tref3)
            call apply_truncation(xa)
            call apply_truncation(xb)
            call apply_truncation(xc)
            call apply_truncation(xd)
            call apply_truncation(xe)
            call apply_truncation(xf)
            call apply_truncation(xg)
            call apply_truncation(xh)
            call apply_truncation(xj)
            call apply_truncation(dhsx)
            call apply_truncation(elz)
        end subroutine truncate_dyncon2
end module
