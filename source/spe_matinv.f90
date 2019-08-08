subroutine ludcmp(a,n,np,indx,d)
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: n, np
    real(dp), intent(inout) :: a(np,np), d
    integer, intent(inout) :: indx(n)

    integer, parameter :: nmax = 100
    real(dp), parameter :: tiniest = 1.0e-20
    integer :: i, j, k, imax
    real(dp) :: vv(nmax), aamax, dum, accum

    d = 1.0_dp

    do i=1,n
        aamax=0.0_dp
        do j=1,n
            if(abs(a(i,j))>aamax) aamax=abs(a(i,j))
        end do
        if(aamax==0.0_dp) stop 'singular'
        vv(i)=1.0_dp/aamax
    end do

    do j=1,n
        if(j>1) then
            do i=1,j-1
                accum=a(i,j)
                if(i>1) then
                    do k=1,i-1
                        accum=accum-a(i,k)*a(k,j)
                    end do
                    a(i,j)=accum
                end if
            end do
        end if

        aamax=0.0_dp
        do i=j,n
            accum=a(i,j)
            if(j>1) then
                do k=1,j-1
                  accum=accum-a(i,k)*a(k,j)
                end do
                a(i,j)=accum
            end if
            dum=vv(i)*abs(accum)
            if(dum>=aamax) then
                imax=i
                aamax=dum
            end if
        end do

        if(j/=imax) then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
        end if

        indx(j)=imax
        if(j/=n) then
            if(a(j,j)==0) a(j,j)=tiniest
            dum=1.0_dp/a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            end do
        end if
    end do

    if(a(n,n)==0.0_dp) a(n,n)=tiniest
end subroutine ludcmp

subroutine lubksb(a,n,np,indx,b)
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: n, np, indx(n)
    real(dp), intent(inout) :: a(np,np), b(n)

    integer :: ii, i, ll, j
    real(dp) :: accum

    ii=0

    do i=1,n
        ll=indx(i)
        accum=b(ll)
        b(ll)=b(i)
        if(ii/=0) then
            do j=ii,i-1
                accum=accum-a(i,j)*b(j)
            end do
        else if(accum /=0) then
            ii=i
        end if
        b(i)=accum
    end do

    do i=n,1,-1
        accum=b(i)
        if(i<n) then
          do j=i+1,n
            accum=accum-a(i,j)*b(j)
          end do
        end if
        b(i)=accum/a(i,i)
    end do
end subroutine lubksb

subroutine inv(a,y,indx,n)
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: n
    real(dp), intent(inout) :: a(n,n), y(n,n)
    integer, intent(inout) :: indx(n)

    integer :: i
    real(dp) :: d

    y = 0.0_dp

    do i=1,n
        y(i,i)=1.0_dp
    end do

    call ludcmp(a,n,n,indx,d)

    do i=1,n
        call lubksb(a,n,n,indx,y(1,i))
    end do
end subroutine inv
