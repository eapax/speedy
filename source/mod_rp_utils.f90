module mod_rp_utils
    use rp_emulator

    implicit none

    private
    public upack1, upack2

    contains
        function upack1(x) result(z)
            type(rpe_complex_var), intent(in) :: x(:,:)
            type(rpe_var) :: z(2*size(x,1), size(x,2))
            integer :: i, j

            do i = 1, size(x,2)
                do j = 1, size(x,1)
                    z(2*j-1,i) = realpart(x(j,i))
                    z(2*j,i)   = imagpart(x(j,i))
                end do
            end do
        end function upack1

        function upack2(x) result(z)
            type(rpe_complex_var), intent(in) :: x(:,:)
            type(rpe_var) :: z(2, size(x,1), size(x,2))
            integer :: i, j

            do i = 1, size(x,2)
                do j = 1, size(x,1)
                    z(1,j,i) = realpart(x(j,i))
                    z(2,j,i) = imagpart(x(j,i))
                end do
            end do
        end function upack2
end module mod_rp_utils
