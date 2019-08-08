module mod_solar

    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    private
    public setup_solar_forcing, sol_oz
    public fsol, ozone, ozupp, zenit, stratz

    namelist /solar_forcing/ solc, epssw

    ! Radiation and cloud constants
    ! solc   = Solar constant (area averaged) in W/m^2
    real(dp) :: solc

    ! epssw  = fraction of incoming solar radiation absorbed by ozone
    real(dp) :: epssw

    ! Zonally-averaged fields for SW/LW scheme (updated in sol_oz)
    ! fsol   = flux of incoming solar radiation
    ! ozone  = flux absorbed by ozone (lower stratos.)
    ! ozupp  = flux absorbed by ozone (upper stratos.)
    ! zenit  = optical depth ratio (function of solar zenith angle)
    ! stratz = stratospheric correction for polar night
    type(rpe_var), dimension(:), allocatable :: &
            fsol, ozone, ozupp, zenit, stratz

    contains
        subroutine setup_solar_forcing(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, solar_forcing)
            write(*, solar_forcing)

            allocate(fsol(ngp))
            allocate(ozone(ngp))
            allocate(ozupp(ngp))
            allocate(zenit(ngp))
            allocate(stratz(ngp))
        end subroutine setup_solar_forcing

        subroutine sol_oz(tyear)
            !  subroutine sol_oz (tyear)
            !
            !  Purpose: Compute zonally-averaged fields to be used
            !           in the computation of SW absorption:
            !           fsol   = flux of incoming solar radiation
            !           ozone  = flux absorbed by ozone (lower stratos.)
            !           ozupp  = flux absorbed by ozone (upper stratos.)
            !           zenit  = function of solar zenith angle
            !  Input:   tyear  = time as fraction of year (0-1, 0 = 1jan.h00)
            !  Updated common blocks: radzon
            use mod_physcon, only: slat, clat

            real(dp), intent(in) :: tyear
            real(dp) :: topsr(il), alpha, azen, coz1, coz2, czen, &
                    dalpha, flat2, fs0
            real(dp) :: nzen, rzen, szen
            integer :: i, j, j0

            ! alpha = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
            alpha=4.0_dp*asin(1.0_dp)*&
                    (tyear+10.0_dp/365.0_dp)
            dalpha=0.0_dp

            coz1= 1.0_dp*max(0.0_dp,cos(alpha-dalpha))
            coz2= 1.8_dp

            azen=1.0_dp
            nzen=2

            rzen=-cos(alpha)*23.45_dp*&
                    asin(1.0_dp)/90.0_dp
            czen=cos(rzen)
            szen=sin(rzen)

            fs0=6.0_dp

            ! Solar radiation at the top
            call solar(tyear,4.0_dp*solc,il,clat,slat,topsr)

            do j=1,il
                j0=1+ix*(j-1)
                flat2=1.5_dp*slat(j)**2-0.5_dp

                ! Solar radiation at the top
                fsol(j0)=topsr(j)

                ! Ozone depth in upper and lower stratosphere
                ozupp(j0)=0.5_dp*epssw
                ozone(j0)=0.4_dp*epssw*&
                        (1.0_dp+coz1*slat(j)+coz2*flat2)

                ! Zenith angle correction to (downward) absorptivity
                zenit(j0)=1.0_dp + &
                        azen*(1.0_dp - &
                                (clat(j)*czen+slat(j)*szen))**nzen

                ! Ozone absorption in upper and lower stratosphere
                ozupp(j0)=fsol(j0)*ozupp(j0)*zenit(j0)
                ozone(j0)=fsol(j0)*ozone(j0)*zenit(j0)

                ! Polar night cooling in the stratosphere
                stratz(j0)=max(fs0-fsol(j0),0.0_dp)

                do i=1,ix-1
                    fsol  (i+j0) = fsol  (j0)
                    ozone (i+j0) = ozone (j0)
                    ozupp (i+j0) = ozupp (j0)
                    zenit (i+j0) = zenit (j0)
                    stratz(i+j0) = stratz(j0)
                end do
            end do

            ! Truncate output for radsw
            call set_precision('Short-Wave Radiation')
            call apply_truncation(fsol)
            call apply_truncation(ozone)
            call apply_truncation(ozupp)
            call apply_truncation(zenit)
            call apply_truncation(stratz)
            call set_precision('Double')
        end subroutine sol_oz

        subroutine solar(tyear,csol,il,clat,slat,topsr)
            ! Average daily flux of solar radiation, from Hartmann (1994)

            real(dp), intent(in) :: tyear, csol
            integer, intent(in) :: il
            real(dp), dimension(il), intent(in) :: clat, slat
            real(dp), intent(inout) :: topsr(il)

            integer :: j
            real(dp) :: ca1, ca2, ca3, cdecl, ch0, csolp, decl, fdis, &
                    h0, alpha, pigr, sa1
            real(dp) :: sa2, sa3, sdecl, sh0, tdecl

            ! 1. Compute declination angle and Earth-Sun distance factor
            pigr  = 2.0_dp*asin(1.0_dp)
            alpha = 2.0_dp*pigr*tyear

            ca1 = cos(alpha)
            sa1 = sin(alpha)
            ca2 = ca1*ca1-sa1*sa1
            sa2 = 2.0_dp*sa1*ca1
            ca3 = ca1*ca2-sa1*sa2
            sa3 = sa1*ca2+sa2*ca1

            decl = 0.006918_dp - &
                    0.399912_dp*ca1 + &
                    0.070257_dp*sa1 - &
                    0.006758_dp*ca2 + &
                    0.000907_dp*sa2 - &
                    0.002697_dp*ca3 + &
                    0.001480_dp*sa3

            fdis = 1.000110_dp + &
                    0.034221_dp*ca1 + &
                    0.001280_dp*sa1 + &
                    0.000719_dp*ca2 + &
                    0.000077_dp*sa2

            cdecl = cos(decl)
            sdecl = sin(decl)
            tdecl = sdecl/cdecl

            ! 2. Compute daily-average insolation at the atm. top
            csolp=csol/pigr

            do j=1,il
                ch0 = min(1.0_dp, &
                        max(-1.0_dp,-tdecl*slat(j)/clat(j)))
                h0  = acos(ch0)
                sh0 = sin(h0)

                topsr(j) = csolp*fdis*(h0*slat(j)*sdecl+sh0*clat(j)*cdecl)
            end do
        end subroutine solar

end module mod_solar
