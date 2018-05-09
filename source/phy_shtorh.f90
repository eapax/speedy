module humidity

    use rp_emulator
    use mod_prec

    implicit none

    private
    public shtorh, q_sat, shtorh_dp

    contains

        subroutine shtorh(imode,ngp,ta,ps,sig,qa,rh,qsat)
            ! subroutine shtorh (imode,ngp,ta,ps,sig,qa,rh,qsat)
            !
            ! Purpose: compute saturation specific humidity and
            !          relative hum. from specific hum. (or viceversa)
            ! Input:  imode : mode of operation
            !         ngp   : no. of grid-points
            !         ta    : abs. temperature
            !         ps    : normalized pressure   (=  p/1000_hPa) [if sig < 0]
            !               : normalized sfc. pres. (= ps/1000_hPa) [if sig > 0]
            !         sig   : sigma level
            !         qa    : specific humidity in g/kg [if imode > 0]
            !         rh    : relative humidity         [if imode < 0]
            !         qsat  : saturation spec. hum. in g/kg
            ! Output: rh    : relative humidity         [if imode > 0]
            !         qa    : specific humidity in g/kg [if imode < 0]

            integer, intent(in) :: imode, ngp
            type(rpe_var), intent(in) :: ta(ngp), ps(*), sig
            type(rpe_var) :: qsat(ngp), qa(ngp), rh(ngp)

            ! 2. Compute rel.hum. RH=Q/Qsat (imode>0), or Q=RH*Qsat (imode<0)
            qsat = q_sat(ngp, ta, ps, sig)
            if (imode.gt.0) then
                rh=qa/qsat
            else if (imode.lt.0) then
                qa=rh*qsat
            end if
        end subroutine shtorh

        function q_sat(ngp, ta, ps, sig)
            ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres.
            ! P (= p/1000_hPa)
            ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.

            integer, intent(in) :: ngp
            type(rpe_var), intent(in) :: ta(ngp), ps(*), sig
            type(rpe_var) :: e0, c1, c2, t0, t1, t2
            type(rpe_var) :: q_sat(ngp)

            integer :: j

            e0 = 6.108e-3
            c1 = 17.269
            c2 = 21.875
            t0 = 273.16
            t1 = 35.86
            t2 = 7.66

            do j=1,ngp
                if (ta(j).ge.t0) then
                    ! Saturation relative to liquid water
                    q_sat(j)=e0*exp(c1*(ta(j)-t0)/(ta(j)-t1))
                else
                    ! Saturation relative to ice
                    q_sat(j)=e0*exp(c2*(ta(j)-t0)/(ta(j)-t2))
                end if
            end do

            if (sig.le.0.0) then
                do j=1,ngp
                    q_sat(j)=rpe_literal(622.)*q_sat(j)/(ps(1)-rpe_literal(0.378)*q_sat(j))
                end do
            else
                do j=1,ngp
                    q_sat(j)=rpe_literal(622.)*q_sat(j)/(sig*ps(j)-rpe_literal(0.378)*q_sat(j))
                end do
            end if

        end function q_sat

        subroutine shtorh_dp(imode,ngp,ta,ps,sig,qa,rh,qsat)
            ! subroutine shtorh (imode,ngp,ta,ps,sig,qa,rh,qsat)
            !
            ! Purpose: compute saturation specific humidity and
            !          relative hum. from specific hum. (or viceversa)
            ! Input:  imode : mode of operation
            !         ngp   : no. of grid-points
            !         ta    : abs. temperature
            !         ps    : normalized pressure   (=  p/1000_hPa) [if sig < 0]
            !               : normalized sfc. pres. (= ps/1000_hPa) [if sig > 0]
            !         sig   : sigma level
            !         qa    : specific humidity in g/kg [if imode > 0]
            !         rh    : relative humidity         [if imode < 0]
            !         qsat  : saturation spec. hum. in g/kg
            ! Output: rh    : relative humidity         [if imode > 0]
            !         qa    : specific humidity in g/kg [if imode < 0]

            integer, intent(in) :: imode, ngp
            real(dp), intent(in) :: ta(ngp), ps(*), sig
            real(dp) :: qsat(ngp)
            real(dp) :: qa(ngp), rh(ngp)

            ! 2. Compute rel.hum. RH=Q/Qsat (imode>0), or Q=RH*Qsat (imode<0)
            qsat = q_sat_dp(ngp, ta, ps, sig)
            if (imode.gt.0) then
                rh=qa/qsat
            else if (imode.lt.0) then
                qa=rh*qsat
            end if
        end subroutine shtorh_dp

        function q_sat_dp(ngp, ta, ps, sig)
            ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres.
            ! P (= p/1000_hPa)
            ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.

            integer, intent(in) :: ngp
            real(dp), intent(in) :: ta(ngp), ps(*), sig
            real(dp) :: e0, c1, c2, t0, t1, t2
            real(dp) :: q_sat_dp(ngp)

            integer :: j

            e0 = 6.108e-3
            c1 = 17.269
            c2 = 21.875
            t0 = 273.16
            t1 = 35.86
            t2 = 7.66

            do j=1,ngp
                if (ta(j).ge.t0) then
                  q_sat_dp(j)=e0*exp(c1*(ta(j)-t0)/(ta(j)-t1))
                else
                  q_sat_dp(j)=e0*exp(c2*(ta(j)-t0)/(ta(j)-t2))
                end if
            end do

            if (sig.le.0.0) then
                do j=1,ngp
                    q_sat_dp(j)=622.*q_sat_dp(j)/(ps(1)-0.378*q_sat_dp(j))
                end do
            else
                do j=1,ngp
                    q_sat_dp(j)=622.*q_sat_dp(j)/(sig*ps(j)-0.378*q_sat_dp(j))
                end do
            end if

        end function q_sat_dp

        subroutine zmeddy(nlon,nlat,ff,zm,eddy)
            ! Decompose a field into zonal-mean and eddy component

            integer, intent(in) :: nlon, nlat
            type(rpe_var), intent(in) :: ff(nlon,nlat)
            type(rpe_var), intent(inout) :: zm(nlat), eddy(nlon,nlat)
            integer :: i, j
            type(rpe_var) :: rnlon

            rnlon=1./nlon

            do j=1,nlat
                zm(j)=0.
                do i=1,nlon
                    zm(j)=zm(j)+ff(i,j)
                end do
                zm(j)=zm(j)*rnlon

                do i=1,nlon
                    eddy(i,j)=ff(i,j)-zm(j)
                end do
            end do
        end subroutine zmeddy
end module humidity
