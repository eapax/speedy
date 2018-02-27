module mod_atparam
    implicit none

    namelist /resolution/ ntrun, mtrun, ix, iy, kx

    integer :: isc = 1
    integer :: ntrun = 30, mtrun = 30, ix = 96, iy = 24, kx=8
    integer :: nx, mx, mxnx, mx2
    integer :: il, ngp, ntrun1
    integer :: nxp, mxp, lmax
    integer :: kx2, kxm, kxp, ntr=1

    contains
        subroutine setup_resolution(fid)
            integer, intent(in) :: fid

            read(fid, resolution)

            nx = ntrun+2
            mx = mtrun+1
            mxnx = mx*nx
            mx2 = 2*mx
            il = 2*iy
            ngp = ix*il
            ntrun1 = ntrun+1
            nxp = nx+1
            mxp = isc*mtrun+1
            lmax = mxp+nx-2
            kx2=2*kx
            kxm=kx-1
            kxp=kx+1
            ntr=1

            write(*, resolution)
        end subroutine setup_resolution
end module
