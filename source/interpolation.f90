module interpolation

    implicit none

    private
    public linear_interp_weights, linear_interp

    contains

        ! Linear interpolation of 1d field: y = f(x)
        ! y(x+dx) = y(x) + dx * dy/dx
        ! dy/dx ~ (y(i) - y(i-1)) / (x(i) - x(i-1))
        ! y(x+dx) ~ a*y(i) + (1-a)*y(i-1)
        ! where a = dx / (x(i) - x(i-1))

        ! Calculate weightings for linear interpolation from x_in to x_out
        subroutine linear_interp_weights( &
                x_in, x_out, nx_in, nx_out, period, idx, weights)

            ! 1d monotononically increasing coordinates of input and output
            integer, intent(in) :: nx_in, nx_out
            real, intent(in) :: x_in(nx_in), x_out(nx_out)

            ! Period for wrapping coordinate (e.g. longitude: period=360)
            real, intent(in) :: period

            ! One index for each element in x_out with the fractional index to
            ! interpolate to in x_in
            integer, intent(out) :: idx(nx_out)
            real, intent(out) :: weights(nx_out)

            integer :: i, n

            i = 1
            do n=1,nx_out
                ! Find the first index where x_out > x_in. i.e. x_out(n) is
                ! between x_in(i-1) and x_in(i)
                do while (x_in(i) < x_out(n) .and. i <= nx_in)
                    i = i+1
                end do

                ! Handle boundary issues
                if (i==1) then
                    if (period > 0) then
                        idx(n) = 1
                        weights(n) = (x_out(n) - x_in(nx_in) + period) / &
                                     (x_in (1) - x_in(nx_in) + period)
                    else
                        ! Assume continous outside boundary: y(1-...) = y(1)
                        idx(n) = 2
                        weights(n) = 0.
                    end if
                else if (i > nx_in) then
                    if (period > 0) then
                        idx(n) = 1
                        weights(n) = (x_out(n) - x_in(nx_in) + period) / &
                                     (x_in (1) - x_in(nx_in) + period)
                    else
                        ! Assume continous outside boundary: y(i+...) = y(i)
                        idx(n) = nx_in
                        weights(n) = 1.
                    end if
                else
                    idx(n) = i
                    weights(n) = (x_out(n) - x_in(i-1)) / (x_in (i) - x_in(i-1))
                end if
            end do

        end subroutine linear_interp_weights

        ! Horizontal interpolation between two lon/lat grids
        function linear_interp(x_in, idx, weights, nx_in, nx_out) result(x_out)
            integer, intent(in) :: nx_in, nx_out
            integer, intent(in) :: idx(nx_out)
            real, intent(in) :: weights(nx_out)
            real, intent(in) :: x_in(nx_in)
            real :: x_out(nx_out)
            integer :: idxp1, idxm1
            integer :: n

            do n=1, nx_out
                ! Check for wrapping coordinates in index
                if (idx(n) == 1) then
                    idxp1 = 1
                    idxm1 = nx_in
                else
                    idxp1 = idx(n)
                    idxm1 = idx(n) - 1
                end if
                x_out(n) =    weights(n)  * x_in(idxp1)   + &
                           (1-weights(n)) * x_in(idxm1)
            end do
        end function linear_interp
end module interpolation