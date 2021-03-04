module mod_dyncon2
    use mod_atparam
    use rp_emulator

    implicit none

    ! Temp. profile for semi-imp. scheme (initial. in IMPINT)
    type(rpe_var), dimension(:),       allocatable :: tref, tref1, tref2, tref3
    type(rpe_var), dimension(:, :),    allocatable :: xc, xd, xe
    type(rpe_var), dimension(:, :, :), allocatable :: xf, xj
    type(rpe_var), allocatable :: dhsx(:), elz(:, :)

    contains
        subroutine setup_dyncon2()
            allocate(tref(kx))
            allocate(tref1(kx))
            allocate(tref2(kx))
            allocate(tref3(kx))

            allocate(xc(kx, kx))
            allocate(xd(kx, kx))
            allocate(xe(kx, kx))

            allocate(xj(kx, kx, lmax))

            allocate(dhsx(kx))
            allocate(elz(mx,nx))
        end subroutine setup_dyncon2

        subroutine truncate_dyncon2()
            ! Todo - truncate tref variables used in other routines than implic
            call apply_truncation(tref)
!            call apply_truncation(tref1)
!            call apply_truncation(tref2)
!            call apply_truncation(tref3)
!
!            call apply_truncation(xc)
!            call apply_truncation(xd)
!            call apply_truncation(xe)
!            call apply_truncation(xj)
!            call apply_truncation(dhsx)
!            call apply_truncation(elz)
        end subroutine truncate_dyncon2
end module
