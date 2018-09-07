module ppo_plevs

    use mod_prec, only: sp, dp

    implicit none

    private
    public pressure_levels, np

    integer, parameter :: np=8

    contains
        ! Interpolate a variable from sigma levels to pressure levels
        subroutine pressure_levels(varid, x_sigma, x_pressure)
            use mod_atparam
            use mod_dynvar, only: ps, T
            use mod_physvar, only: tg1
            use mod_physcon, only: sigl, rd, gg, pout
            use spectral, only: grid

            integer, intent(in) :: varid
            real(dp), dimension(ngp, kx), intent(in) :: x_sigma
            real(sp), dimension(ngp, kx), intent(out) :: x_pressure

            real(dp), dimension(ngp) :: x_pressure_dp
            real(dp), dimension(kx) :: zinp, rdzinp
            integer :: k0(ngp)
            real(dp) :: w0(ngp)
            real(dp), dimension(ngp) :: psgr
            real(dp), dimension(ngp) :: zout
            real(dp), dimension(ngp) :: T_pressure
            real(dp) :: textr, aref, tref, phi1, phi2
            integer :: k, j

            ! Vertical interpolation from sigma level to pressure level
            ! sigl is constant so this should only be done once
            zinp = -sigl
            do k=2,kx
               rdzinp(k) = 1.0_dp/(zinp(k-1)-zinp(k))
            end do

            ! This should only be done once. Not for each variable
            call grid(ps, psgr, 1)

            do k=1,np
                ! psgr changes each timestep so this should always be recalculated
                ! This should only be done once. Not for each variable
                do j=1,ngp
                    zout(j) = psgr(j) - log(pout(k))
                end do

                call setvin(zinp,rdzinp,zout,ngp,kx,k0,w0)

                ! This is done for all variables apart from temperature
                if (varid/=3 .and. varid/=5) then
                    do j=1,ngp
                        w0(j) = max(w0(j),0.0_dp)
                    end do
                end if

                ! Do the interpolation
                call verint(x_pressure_dp,x_sigma,ngp,kx,k0,w0)

                ! Pass variable to single precision output
                if (varid==3) then
                    ! Corrections applied to temperature
                    do j=1,ngp
                        if(zout(j)<zinp(kx)) then
                            textr = max(x_pressure_dp(j), x_sigma(j,kx))
                            aref = rd*0.006_dp/gg * (zinp(kx)-zout(j))
                            tref = x_sigma(j,kx)*(1.0_dp+aref+0.5_dp*aref*aref)
                            x_pressure_dp(j) = textr + 0.7_dp*(tref-textr)
                        end if
                    end do
                    x_pressure(:,k) = x_pressure_dp

                else if (varid==5) then
                    call verint(T_pressure,Tg1,ngp,kx,k0,w0)
                    ! Corrections applied to temperature
                    do j=1,ngp
                        if(zout(j)<zinp(kx)) then
                            textr = max(T_pressure(j), Tg1(j,kx))
                            aref = rd*0.006_dp/gg * (zinp(kx)-zout(j))
                            tref = Tg1(j,kx)*(1.0_dp+aref+0.5_dp*aref*aref)
                            T_pressure(j) = textr + 0.7_dp*(tref-textr)
                        end if
                    end do

                    do j=1,ngp
                        w0(j) = max(w0(j),0.0_dp)
                    end do

                    ! Corrections applied to geopotential height
                    do j=1,ngp
                        phi1 = x_sigma(j,k0(j)) &
                           & +0.5_dp*rd*(T_pressure(j)+Tg1(j,k0(j)))* &
                                        (zout(j)-zinp(k0(j)))
                        phi2 = x_sigma(j,k0(j)-1) &
                           & +0.5_dp*rd*(T_pressure(j)+Tg1(j,k0(j)-1))* &
                                        (zout(j)-zinp(k0(j)-1))
                        x_pressure_dp(j) = phi1 + w0(j)*(phi2-phi1)
                    end do
                    x_pressure(:,k) = x_pressure_dp / gg
                else
                    x_pressure(:,k) = x_pressure_dp
                end if

            end do
        end subroutine pressure_levels

        subroutine setvin(zinp,rdzinp,zout,ngp,nlev,k0,w0)
            implicit none

            integer, intent(in) :: ngp, nlev
            real(dp), intent(in) :: zinp(nlev), rdzinp(nlev), zout(ngp)
            real(dp), intent(out) :: w0(ngp)
            integer, intent(out) :: k0(ngp)
            integer :: j, k

            ! *** 1. Select closest vertical levels
            do j=1,ngp
                k0(j)=2
            end do

            do k=2,nlev-1
                do j=1,ngp
                    if (zout(j)<zinp(k)) k0(j)=k+1
                end do
            end do

            ! *** 2. Compute interpolation weight
            do j=1,ngp
                w0(j)=(zout(j)-zinp(k0(j)))*rdzinp(k0(j))
            end do
        end subroutine setvin

        subroutine verint(f2d,f3d,ngp,nlev,k0,w0)
            implicit none

            ! *** 1. Perform vertical interpolation
            integer, intent(in) :: ngp, nlev, k0(ngp)
            real(dp), intent(in) :: f3d(ngp,nlev), w0(ngp)
            real(dp), intent(out) :: f2d(ngp)
            integer :: j

            do j=1,ngp
                f2d(j)=f3d(j,k0(j))+w0(j)*(f3d(j,k0(j)-1)-f3d(j,k0(j)))
            end do
        end subroutine verint

end module ppo_plevs