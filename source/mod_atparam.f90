module mod_atparam
    implicit none

    namelist /resolution/ ntrun, mtrun, ix, iy, kx

    integer :: isc = 1

    ! Model resolution (read in from namelist)
    integer :: ntrun = 30, mtrun = 30, ix = 96, iy = 24, kx=8

    ! Constants derived from model resolution at setup
    integer :: nx, mx, mx2
    integer :: il, ngp, ntrun1
    integer :: nxp, mxp, lmax
    integer :: kxm, kxp

    ! Number of tracers (1 for humidity)
    integer :: ntr

    contains
        subroutine setup_resolution(fid)
            integer, intent(in) :: fid

            read(fid, resolution)

            nx = ntrun+2
            mx = mtrun+1
            mx2 = 2*mx
            il = 2*iy
            ngp = ix*il
            ntrun1 = ntrun+1
            nxp = nx+1
            mxp = isc*mtrun+1
            lmax = mxp+nx-2
            kxm=kx-1
            kxp=kx+1
            ntr=1

            write(*, resolution)
        end subroutine setup_resolution
end module
