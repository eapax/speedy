module mod_fordate
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public setup_forcing, ini_fordate, fordate
    public ablco2, albsea, albice
    public alb_l, alb_s, albsfc, snowc

    namelist /forcing/ lco2, iyear_ref, del_co2, ablco2, &
            albsea, albice, albsn

    ! Flag and parameters for co2 forcing trend
    ! Flag for CO2 optical thickness increase
    logical :: lco2
    ! If lco2=.true., the year corresponding to ablco2 baseline
    integer :: iyear_ref
    ! If lco2=.true., the yearly increase in co2 absorbtion from iyear_ref
    real(dp) :: del_co2
    ! ablco2 = abs. of air in CO2 band
    real(dp) :: ablco2, ablco2_ref

    ! Constants for surface albedos
    ! albsea = Albedo over sea
    real(dp) :: albsea
    ! albice = Albedo over sea ice (for ice fraction = 1)
    real(dp) :: albice
    ! albsn  = Albedo over snow (for snow cover = 1)
    real(dp) :: albsn

    ! Radiative properties of the surface (updated in fordate)
    ! alb_l, alb_s and snowc used in phy_suflux. albsfc used in phy_radsw
    ! alb_l  = daily-mean albedo over land (bare-land + snow)
    ! alb_s  = daily-mean albedo over sea  (open sea + sea ice)
    ! albsfc = combined surface albedo (land + sea)
    ! snowc  = effective snow cover (fraction)
    real(dp), dimension(:), allocatable :: alb_l, alb_s, albsfc, snowc

    real(dp) :: gamlat, pexp

    contains
        subroutine setup_forcing(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, forcing)
            write(*, forcing)

            allocate(alb_l(ngp))
            allocate(alb_s(ngp))
            allocate(albsfc(ngp))
            allocate(snowc(ngp))
        end subroutine setup_forcing

        subroutine ini_fordate()
            use mod_dyncon0, only: gamma
            use mod_physcon, only: gg, rd

            ! Calculate local variables
            ablco2_ref = ablco2
            gamlat = gamma/(1000.0_dp * gg)
            pexp = 1.0_dp/(rd * gamlat)

            call fordate()
        end subroutine ini_fordate

        subroutine fordate()
            !
            !   subroutine fordate (imode)
            !
            !   purpose : compute forcing fields for the current date
            !             and correction terms for horiz. diffusion

            use mod_dyncon0, only: refrh1
            use mod_hdifcon, only: tcorh, qcorh
            use mod_surfcon, only: phis0, alb0, sd2sc
            use mod_cli_land, only: fmask_l
            use mod_date, only: iyear, tyear
            use mod_var_land, only: stl_am, snowd_am
            use mod_cli_sea, only: fmask_s
            use mod_var_sea, only: sst_am, sice_am
            use mod_solar, only: sol_oz
            use humidity, only: q_sat
            use spectral, only: spec

            real(dp), dimension(ix, il) :: corh, tsfc, tref, psfc
            real(dp), dimension(ngp) :: qsfc, qref
            real(dp) :: fland(ngp), alb_0(ngp)

            integer :: i, j, ij

            fland = reshape(fmask_l, (/ngp/))
            alb_0 = reshape(alb0, (/ngp/))

            ! time variables for interpolation are set by newdate

            ! 1. daily-mean radiative forcing
            ! incoming solar radiation
            call sol_oz(tyear)

            ! total surface albedo
            do j = 1, ngp
                snowc(j)  = min(1.0_dp, snowd_am(j)/sd2sc)
                alb_l(j)  = alb_0(j) + snowc(j) * (albsn - alb_0(j))
                alb_s(j)  = albsea + sice_am(j) * (albice - albsea)
                albsfc(j) = alb_s(j) + fland(j) * (alb_l(j) - alb_s(j))
            end do

            ! linear trend of co2 absorptivity (del_co2: rate of change per year)
            if (lco2) then
                ablco2 = ablco2_ref * exp(del_co2 * (iyear + tyear - iyear_ref))
            end if

            ! 2. temperature correction term for horizontal diffusion
            corh = gamlat * phis0

            call spec(corh,tcorh)

            ! 3. humidity correction term for horizontal diffusion
            ij = 0
            do j = 1, il
                do i = 1, ix
                    ij = ij + 1
                    tsfc(i,j) = fmask_l(i,j) * stl_am(ij)&
                        & + fmask_s(i,j) * sst_am(ij)
                    tref(i,j) = tsfc(i,j) + corh(i,j)
                    psfc(i,j) = (tsfc(i,j)/tref(i,j))**pexp
                end do
            end do

            qref = q_sat(ngp, reshape(tref, (/ngp/)), &
                    (/1.0_dp/), -1.0_dp)
            qsfc = q_sat(ngp, reshape(tsfc, (/ngp/)), psfc,  1.0_dp)

            corh = refrh1 * reshape((qref - qsfc), (/ix, il/))

            call spec(corh,qcorh)
        end subroutine fordate
end module mod_fordate
