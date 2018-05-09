! Functions and parameters required to start a model run from a different
! resolution restart file
module downscaling
    use mod_atparam, only: ix, il, iy
    use spectral, only: sia, gaussl
    use interpolation, only: linear_interp_weights, linear_interp
    use rp_emulator

    implicit none

    private
    public nx_in, mx_in, ix_in, il_in, kx_in, &
            setup_downscaling, calc_grid_weights, regrid

    namelist /input_resolution/ ntrun_in, mtrun_in, ix_in, iy_in, kx_in

    integer :: ntrun_in, mtrun_in
    integer :: ix_in, iy_in, kx_in
    integer :: nx_in, mx_in, il_in

    integer, allocatable :: idx_x(:), idx_y(:)
    type(rpe_var), allocatable :: weights_x(:), weights_y(:)

    contains
        subroutine setup_downscaling(fid)
            integer, intent(in) :: fid

            read(fid, input_resolution)

            nx_in = ntrun_in+2
            mx_in = mtrun_in+1
            il_in = 2*iy_in

            allocate(idx_x(ix))
            allocate(idx_y(il))
            allocate(weights_x(ix))
            allocate(weights_y(il))

        end subroutine setup_downscaling

        subroutine calc_grid_weights()
            type(rpe_var) :: lons_in(ix_in), lons_out(ix)
            type(rpe_var) :: sin_lat(iy_in), wt_in(iy_in)
            type(rpe_var) :: sin_lat_in(il_in), sin_lat_out(il)
            integer :: n, j, jj

            ! Calculate interpolation weights for longitude
            lons_in =  (/ (n*(360./ix_in), n=1,ix_in) /)
            lons_out = (/ (n*(360./ix)   , n=1,ix   ) /)
            call linear_interp_weights( &
                    lons_in, lons_out, ix_in, ix, rpe_literal(360.), idx_x, weights_x)

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
                    sin_lat_in, sin_lat_out, il_in, il, rpe_literal(0.), idx_y, weights_y)
        end subroutine

        ! Regrid an input field to the model resolution
        subroutine regrid(xgrid_in, xgrid_out)
            ! Input field
            type(rpe_var), intent(in) :: xgrid_in(ix_in, il_in)
            ! Field interpolated to model resolution
            type(rpe_var) :: xgrid_out(ix, il)
            ! Intermediate step (longitude->latitude)
            type(rpe_var) :: xgrid_inter(ix, il_in)
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