module humidity

    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public shtorh, q_sat, shtorh_celsius, q_sat_celsius, zero_c

    ! Conversion constant between Kelvin and Celsius
    real(dp) :: zero_C = 273.0_dp

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
            if (imode>0) then
                rh=qa/qsat
            else if (imode<0) then
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

            e0 = rpe_literal(6.108d-3)
            c1 = rpe_literal(17.269_dp)
            c2 = rpe_literal(21.875_dp)
            t0 = rpe_literal(273.16_dp)
            t1 = rpe_literal(35.86_dp)
            t2 = rpe_literal(7.66_dp)

            do j=1,ngp
                if (ta(j)>=t0) then
                    ! Saturation relative to liquid water
                    q_sat(j)=e0*exp(c1*(ta(j)-t0)/(ta(j)-t1))
                else
                    ! Saturation relative to ice
                    q_sat(j)=e0*exp(c2*(ta(j)-t0)/(ta(j)-t2))
                end if
            end do

            if (sig<=rpe_literal(0.0_dp)) then
                do j=1,ngp
                    q_sat(j)=rpe_literal(622.0_dp)*q_sat(j)/ &
                            (ps(1)-rpe_literal(0.378_dp)*q_sat(j))
                end do
            else
                do j=1,ngp
                    q_sat(j)=rpe_literal(622.0_dp)*q_sat(j)/ &
                            (sig*ps(j)-rpe_literal(0.378_dp)*q_sat(j))
                end do
            end if

        end function q_sat

        subroutine shtorh_celsius(imode,ngp,ta,ps,sig,qa,rh,qsat)
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
            qsat = q_sat_celsius(ngp, ta, ps, sig)
            if (imode>0) then
                rh=qa/qsat
            else if (imode<0) then
                qa=rh*qsat
            end if
        end subroutine shtorh_celsius

        function q_sat_celsius(ngp, ta, ps, sig)
            ! 1. Compute Qsat (g/kg) from T (degK) and normalized pres.
            ! P (= p/1000_hPa)
            ! If sig > 0, P = Ps * sigma, otherwise P = Ps(1) = const.

            integer, intent(in) :: ngp
            type(rpe_var), intent(in) :: ta(ngp), ps(*), sig
            type(rpe_var) :: e0, c1, c2, t0, t1, t2
            type(rpe_var) :: q_sat_celsius(ngp)

            integer :: j

            e0 = rpe_literal(6.108d-3)
            c1 = rpe_literal(17.269_dp)
            c2 = rpe_literal(21.875_dp)
            t0 = rpe_literal(273.16_dp - zero_c)
            t1 = rpe_literal(35.86_dp - zero_c)
            t2 = rpe_literal(7.66_dp - zero_c)

            do j=1,ngp
                if (ta(j)>=t0) then
                    ! Saturation relative to liquid water
                    q_sat_celsius(j)=e0*exp(c1*(ta(j)-t0)/(ta(j)-t1))
                else
                    ! Saturation relative to ice
                    q_sat_celsius(j)=e0*exp(c2*(ta(j)-t0)/(ta(j)-t2))
                end if
            end do

            if (sig<=rpe_literal(0.0_dp)) then
                do j=1,ngp
                    q_sat_celsius(j)=rpe_literal(622.0_dp)*q_sat_celsius(j)/ &
                            (ps(1)-rpe_literal(0.378_dp)*q_sat_celsius(j))
                end do
            else
                do j=1,ngp
                    q_sat_celsius(j)=rpe_literal(622.0_dp)*q_sat_celsius(j)/ &
                            (sig*ps(j)-rpe_literal(0.378_dp)*q_sat_celsius(j))
                end do
            end if

        end function q_sat_celsius
end module humidity
