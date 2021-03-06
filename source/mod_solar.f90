module mod_solar

    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    private
    public setup_solar_forcing, truncate_solar, sol_oz
    public fsol, ozone, ozupp, zenit, stratz

    namelist /solar_forcing/ solc, epssw

    ! Radiation and cloud constants
    ! solc   = Solar constant (area averaged) in W/m^2
    type(rpe_var) :: solc

    ! epssw  = fraction of incoming solar radiation absorbed by ozone
    type(rpe_var) :: epssw

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

        subroutine truncate_solar()
            ! Truncate local variables for calculating solar forcing
            ! Namelist variables
            call apply_truncation(solc)
            call apply_truncation(epssw)

        end subroutine truncate_solar

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

            type(rpe_var), intent(in) :: tyear
            type(rpe_var) :: topsr(il), alpha, azen, coz1, coz2, czen, &
                    dalpha, flat2, fs0
            type(rpe_var) :: nzen, rzen, szen
            integer :: i, j, j0

            ! alpha = year phase ( 0 - 2pi, 0 = winter solstice = 22dec.h00 )
            alpha=rpe_literal(4.0_dp)*asin(rpe_literal(1.0_dp))*&
                    (tyear+rpe_literal(10.0_dp)/rpe_literal(365.0_dp))
            dalpha=0.0_dp

            coz1= rpe_literal(1.0_dp)*max(rpe_literal(0.0_dp),cos(alpha-dalpha))
            coz2= 1.8_dp

            azen=1.0_dp
            nzen=2

            rzen=-cos(alpha)*rpe_literal(23.45_dp)*&
                    asin(rpe_literal(1.0_dp))/rpe_literal(90.0_dp)
            czen=cos(rzen)
            szen=sin(rzen)

            fs0=6.0_dp

            ! Solar radiation at the top
            call solar(tyear,rpe_literal(4.0_dp)*solc,il,clat,slat,topsr)

            do j=1,il
                j0=1+ix*(j-1)
                flat2=rpe_literal(1.5_dp)*slat(j)**2-rpe_literal(0.5_dp)

                ! Solar radiation at the top
                fsol(j0)=topsr(j)

                ! Ozone depth in upper and lower stratosphere
                ozupp(j0)=rpe_literal(0.5_dp)*epssw
                ozone(j0)=rpe_literal(0.4_dp)*epssw*&
                        (rpe_literal(1.0_dp)+coz1*slat(j)+coz2*flat2)

                ! Zenith angle correction to (downward) absorptivity
                zenit(j0)=rpe_literal(1.0_dp) + &
                        azen*(rpe_literal(1.0_dp) - &
                                (clat(j)*czen+slat(j)*szen))**nzen

                ! Ozone absorption in upper and lower stratosphere
                ozupp(j0)=fsol(j0)*ozupp(j0)*zenit(j0)
                ozone(j0)=fsol(j0)*ozone(j0)*zenit(j0)

                ! Polar night cooling in the stratosphere
                stratz(j0)=max(fs0-fsol(j0),rpe_literal(0.0_dp))

                do i=1,ix-1
                    fsol  (i+j0) = fsol  (j0)
                    ozone (i+j0) = ozone (j0)
                    ozupp (i+j0) = ozupp (j0)
                    zenit (i+j0) = zenit (j0)
                    stratz(i+j0) = stratz(j0)
                end do
            end do
        end subroutine sol_oz

        subroutine solar(tyear,csol,il,clat,slat,topsr)
            ! Average daily flux of solar radiation, from Hartmann (1994)

            type(rpe_var), intent(in) :: tyear, csol
            integer, intent(in) :: il
            type(rpe_var), dimension(il), intent(in) :: clat, slat
            type(rpe_var), intent(inout) :: topsr(il)

            integer :: j
            type(rpe_var) :: ca1, ca2, ca3, cdecl, ch0, csolp, decl, fdis, &
                    h0, alpha, pigr, sa1
            type(rpe_var) :: sa2, sa3, sdecl, sh0, tdecl

            ! 1. Compute declination angle and Earth-Sun distance factor
            pigr  = rpe_literal(2.0_dp)*asin(rpe_literal(1.0_dp))
            alpha = rpe_literal(2.0_dp)*pigr*tyear

            ca1 = cos(alpha)
            sa1 = sin(alpha)
            ca2 = ca1*ca1-sa1*sa1
            sa2 = rpe_literal(2.0_dp)*sa1*ca1
            ca3 = ca1*ca2-sa1*sa2
            sa3 = sa1*ca2+sa2*ca1

            decl = rpe_literal(0.006918_dp) - &
                    rpe_literal(0.399912_dp)*ca1 + &
                    rpe_literal(0.070257_dp)*sa1 - &
                    rpe_literal(0.006758_dp)*ca2 + &
                    rpe_literal(0.000907_dp)*sa2 - &
                    rpe_literal(0.002697_dp)*ca3 + &
                    rpe_literal(0.001480_dp)*sa3

            fdis = rpe_literal(1.000110_dp) + &
                    rpe_literal(0.034221_dp)*ca1 + &
                    rpe_literal(0.001280_dp)*sa1 + &
                    rpe_literal(0.000719_dp)*ca2 + &
                    rpe_literal(0.000077_dp)*sa2

            cdecl = cos(decl)
            sdecl = sin(decl)
            tdecl = sdecl/cdecl

            ! 2. Compute daily-average insolation at the atm. top
            csolp=csol/pigr

            do j=1,il
                ch0 = min(rpe_literal(1.0_dp), &
                        max(-rpe_literal(1.0_dp),-tdecl*slat(j)/clat(j)))
                h0  = acos(ch0)
                sh0 = sin(h0)

                topsr(j) = csolp*fdis*(h0*slat(j)*sdecl+sh0*clat(j)*cdecl)
            end do
        end subroutine solar

end module mod_solar
