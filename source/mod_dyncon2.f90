module mod_dyncon2
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Temp. profile for semi-imp. scheme (initial. in IMPINT)
    real(dp), dimension(:),       allocatable :: tref, tref1, tref2, tref3
    real(dp), dimension(:, :),    allocatable :: xa, xb, xc, xd, xe
    real(dp), dimension(:, :, :), allocatable :: xf, xg, xh, xj
    real(dp), allocatable :: dhsx(:), elz(:, :)

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
end module
