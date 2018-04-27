! Functions and parameters required to start a model run from a different
! resolution restart file
module downscaling
    use mod_atparam, only: ix, il, iy
    use spectral, only: sia, gaussl
    use interpolation, only: linear_interp_weights, linear_interp

    implicit none

    private
    public nx_in, mx_in, ix_in, il_in, kx_in, calc_grid_weights, regrid

    integer, parameter :: ntrun_in = 30, mtrun_in = 30
    integer, parameter :: nx_in = ntrun_in+2, mx_in = mtrun_in+1
    integer, parameter :: ix_in = 96, iy_in = 24, kx_in = 8
    integer, parameter :: il_in = 2*iy_in

    integer :: idx_x(ix), idx_y(il)
    real :: weights_x(ix), weights_y(il)

    contains
        subroutine calc_grid_weights()
            real :: lons_in(ix_in), lons_out(ix)
            real :: sin_lat(iy_in), wt_in(iy_in)
            real :: sin_lat_in(il_in), sin_lat_out(il)
            integer :: n, j, jj

            ! Calculate interpolation weights for longitude
            lons_in =  (/ (n*(360./ix_in), n=1,ix_in) /)
            lons_out = (/ (n*(360./ix)   , n=1,ix   ) /)
            call linear_interp_weights( &
                    lons_in, lons_out, ix_in, ix, 360., idx_x, weights_x)

            ! Calculate interpolation weights for sin(latitude)
            call gaussl(sin_lat,wt_in,iy_in)

            ! Convert to full array
            do j=1,iy_in
                jj=il_in+1-j
                sin_lat_in(j) = -sin_lat(j)
                sin_lat_in(jj) = sin_lat(j)
            end do

            do j=1,iy
                jj=il+1-j
                sin_lat_out(j) = -sia(j)
                sin_lat_out(jj) = sia(j)
            end do

            call linear_interp_weights( &
                    sin_lat_in, sin_lat_out, il_in, il, 0., idx_y, weights_y)
        end subroutine

        ! Regrid an input field to the model resolution
        subroutine regrid(xgrid_in, xgrid_out)
            ! Input field
            real, intent(in) :: xgrid_in(ix_in, il_in)
            ! Field interpolated to model resolution
            real :: xgrid_out(ix, il)
            ! Intermediate step (longitude->latitude)
            real :: xgrid_inter(ix, il_in)
            integer :: i, j

            ! Interpolate each longitude band to the new grid
            do j=1, il_in
                xgrid_inter(:, j) = linear_interp(xgrid_in(:, j), &
                        idx_x, weights_x, ix_in, ix)
            end do

            ! Interpolate each latitude band to the new grid
            do i=1, ix
                xgrid_out(i, :) = linear_interp(xgrid_inter(i, :), &
                        idx_y, weights_y, il_in, il)
            end do

        end subroutine regrid
end module downscaling