module phy_convmf
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public convmf, setup_convection, ini_convmf

    ! Variables loaded in by namelist
    namelist /convection/ psmin, trcnv, rhbl, rhil, entmax, smf

    ! Minimum (norm.) sfc. pressure for the occurrence of convection
    real(dp) :: psmin

    ! Time of relaxation (in hours) towards reference state
    real(dp) :: trcnv

    ! Relative hum. threshold in the boundary layer
    real(dp) :: rhbl

    ! Rel. hum. threshold in intermed. layers for secondary mass flux
    real(dp) :: rhil

    ! Max. entrainment as a fraction of cloud-base mass flux
    real(dp) :: entmax

    ! Ratio between secondary and primary mass flux at cloud-base
    real(dp) :: smf

    ! Local derived variables
    real(dp) :: fm0
    real(dp), allocatable :: entr(:)

    contains
        subroutine setup_convection(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, convection)
            write(*, convection)

            allocate(entr(2:kx-1))
        end subroutine setup_convection

        subroutine ini_convmf()
            ! Calculate local variables for convection scheme
            use mod_physcon, only: p0, gg, sig, dsig

            integer :: k
            real(dp) :: sentr

            fm0=p0*dsig(kx)/(gg*trcnv*3600.0_dp)

            ! Entrainment profile (up to sigma = 0.5)
            sentr=0.0_dp
            do k=2,kxm
                entr(k)=(max(0.0_dp, sig(k)-0.5_dp))**2
                sentr=sentr+entr(k)
            end do

            sentr=entmax/sentr
            entr(2:kxm) = entr(2:kxm) * sentr
        end subroutine ini_convmf

        subroutine convmf (&
                psa, se, qa, qsat, hflx2tend, flx2tend, &
                itop, cbmf, precnv, dfse, dfqa)
            ! SUBROUTINE CONVMF (PSA,SE,QA,QSAT, ITOP,CBMF,PRECNV,DFSE,DFQA)
            !
            ! Purpose: Compute convective fluxes of dry static energy and
            !          moisture using a simplified mass-flux scheme
            use mod_physcon, only: alhc, sig, wvi

            ! Input:  PSA    = norm. surface pressure [p/p0]            (2-dim)
            real(dp), intent(in) :: psa(ngp)
            !         SE     = dry static energy                        (3-dim)
            real(dp), intent(in) :: se(ngp, kx)
            !         QA     = specific humidity [g/kg]                 (3-dim)
            real(dp), intent(in) :: qa(ngp, kx)
            !         QSAT   = saturation spec. hum. [g/kg]             (3-dim)
            real(dp), intent(in) :: qsat(ngp, kx)
            !         hflx2tend = Conversion factor between heat fluxes and T tendency
            real(dp), intent(in) :: hflx2tend(ngp,kx)
            !         flx2tend = Conversion factor between fluxes and tendencies
            real(dp), intent(in) :: flx2tend(ngp,kx)

            ! Output: ITOP   = top of convection (layer index)          (2-dim)
            integer, intent(out) :: itop(ngp)
            !         CBMF   = cloud-base mass flux                     (2-dim)
            real(dp), intent(out) :: cbmf(ngp)
            !         PRECNV = convective precipitation [g/(m^2 s)]     (2-dim)
            real(dp), intent(out) :: precnv(ngp)
            !         DFSE   = net flux of d.s.en. into each atm. layer (3-dim)
            real(dp), intent(out) :: dfse(ngp,kx)
            !         DFQA   = net flux of sp.hum. into each atm. layer (3-dim)
            real(dp), intent(out) :: dfqa(ngp,kx)

            ! Local variables
            integer :: j, k, k1
            real(dp) :: qdif(ngp), delq, enmass, fdq, fds, fmass, fpsa, fqmax, &
                    fsq, fuq, fus, qb, qmax, qsatb, rdps, sb

            ! 1. Initialization of output and workspace arrays
            fqmax=5.0_dp

            rdps=2.0_dp/(1.0_dp-psmin)

            dfse = 0.0_dp
            dfqa = 0.0_dp

            cbmf = 0.0_dp
            precnv = 0.0_dp

            ! 2. Check of conditions for convection
            call diagnose_convection(psa, se, qa, qsat, itop, qdif)

            ! 3. Convection over selected grid-points
            do j=1,ngp
                if (itop(j)==kxp) cycle

                ! 3.1 Boundary layer (cloud base)
                k =kx
                k1=k-1

                ! Maximum specific humidity in the PBL
                qmax=max(1.01_dp*qa(j,k),qsat(j,k))

                ! Dry static energy and moisture at upper boundary
                sb=se(j,k1)+wvi(k1,2)*(se(j,k)-se(j,k1))
                qb=qa(j,k1)+wvi(k1,2)*(qa(j,k)-qa(j,k1))
                qb=min(qb,qa(j,k))

                ! Cloud-base mass flux, computed to satisfy:
                ! fmass*(qmax-qb)*(g/dp)=qdif/trcnv
                fpsa=psa(j)*min(1.0_dp,(psa(j)-psmin)*rdps)
                fmass=fm0*fpsa*min(fqmax,qdif(j)/(qmax-qb))
                cbmf(j)=fmass

                ! Upward fluxes at upper boundary
                fus=fmass*se(j,k)
                fuq=fmass*qmax

                ! Downward fluxes at upper boundary
                fds=fmass*sb
                fdq=fmass*qb

                ! Net flux of dry static energy and moisture
                dfse(j,k)=fds-fus
                dfqa(j,k)=fdq-fuq

                ! 3.2 Intermediate layers (entrainment)
                do k=kx-1,itop(j)+1,-1
                    k1=k-1

                    ! Fluxes at lower boundary
                    dfse(j,k)=fus-fds
                    dfqa(j,k)=fuq-fdq

                    ! Mass entrainment
                    enmass=entr(k)*psa(j)*cbmf(j)
                    fmass=fmass+enmass

                    ! Upward fluxes at upper boundary
                    fus=fus+enmass*se(j,k)
                    fuq=fuq+enmass*qa(j,k)

                    ! Downward fluxes at upper boundary
                    sb=se(j,k1)+wvi(k1,2)*(se(j,k)-se(j,k1))
                    qb=qa(j,k1)+wvi(k1,2)*(qa(j,k)-qa(j,k1))
                    fds=fmass*sb
                    fdq=fmass*qb

                    ! Net flux of dry static energy and moisture
                    dfse(j,k)=dfse(j,k)+fds-fus
                    dfqa(j,k)=dfqa(j,k)+fdq-fuq

                    ! Secondary moisture flux
                    delq=rhil*qsat(j,k)-qa(j,k)
                    if (delq>0.0_dp) then
                        fsq=smf*cbmf(j)*delq
                        dfqa(j,k)   =dfqa(j,k)   +fsq
                        dfqa(j,kx)=dfqa(j,kx)-fsq
                    end if
                end do

                ! 3.3 Top layer (condensation and detrainment)
                k=itop(j)

                ! Flux of convective precipitation
                qsatb=qsat(j,k)+wvi(k,2)*(qsat(j,k+1)-qsat(j,k))
                precnv(j)=max(fuq-fmass*qsatb,0.0_dp)

                ! Net flux of dry static energy and moisture
                dfse(j,k)=fus-fds+alhc*precnv(j)
                dfqa(j,k)=fuq-fdq-precnv(j)
            end do

            ! Convert fluxes to temperature tendencies
            dfse = dfse*hflx2tend
            dfqa = dfqa*flx2tend
        end subroutine convmf

        subroutine diagnose_convection(psa, se, qa, qsat, itop, qdif)
            ! Purpose: Check criteria for convection in each grid column and
            !          determine the level of convection. Calculate the humidity
            !          excess in convective gridboxes.

            use mod_physcon, only: alhc, wvi

            ! Input:  PSA    = norm. surface pressure [p/p0]            (2-dim)
            real(dp), intent(in) :: psa(ngp)
            !         SE     = dry static energy                        (3-dim)
            real(dp), intent(in) :: se(ngp, kx)
            !         QA     = specific humidity [g/kg]                 (3-dim)
            real(dp), intent(in) :: qa(ngp, kx)
            !         QSAT   = saturation spec. hum. [g/kg]             (3-dim)
            real(dp), intent(in) :: qsat(ngp, kx)

            ! Output: ITOP   = top of convection (layer index)          (2-dim)
            integer, intent(out) :: itop(ngp)
            !         QDIF   = Humidity anomaly in convection gridpoints (2-dim)
            real(dp), intent(out) :: qdif(ngp)

            ! Local variables
            integer :: j, k, ktop1, ktop2
            real(dp) :: mss(ngp,2:kx), mse0, mse1, mss0, mss2, msthr, &
                    qthr0, qthr1, rlhc
            logical :: lqthr

            rlhc=1.0_dp/alhc

            ! Saturation moist static energy
            do k=2,kx
                do j=1,ngp
                    mss(j,k)=se(j,k)+alhc*qsat(j,k)
                end do
            end do

            do j=1,ngp
                itop(j)=kxp

                if (psa(j)>psmin) then
                    ! Minimum of moist static energy in the lowest two levels
                    mse0=se(j,kx)+alhc*qa(j,kx)
                    mse1=se(j,kxm) +alhc*qa(j,kxm)
                    mse1=min(mse0,mse1)

                    ! Saturation (or super-saturated) moist static energy in PBL
                    mss0=max(mse0,mss(j,kx))

                    ktop1=kx
                    ktop2=kx

                    do k=kx-3,3,-1
                        mss2=mss(j,k)+wvi(k,2)*(mss(j,k+1)-mss(j,k))

                        ! Check 1: conditional instability
                        !          (MSS in PBL > MSS at top level)
                        if (mss0>mss2) then
                           ktop1=k
                        end if

                        ! Check 2: gradient of actual moist static energy
                        !          between lower and upper troposphere
                        if (mse1>mss2) then
                           ktop2=k
                           msthr=mss2
                        end if
                    end do

                    if (ktop1<kx) then
                        ! Check 3: RH > RH_c at both k=kx and k=kxm
                        qthr0=rhbl*qsat(j,kx)
                        qthr1=rhbl*qsat(j,kxm)
                        lqthr=(qa(j,kx)>qthr0 .and. qa(j,kxm)>qthr1)

                        if (ktop2<kx) then
                           itop(j)=ktop1
                           qdif(j)=max(qa(j,kx)-qthr0,(mse0-msthr)*rlhc)
                        else if (lqthr) then
                           itop(j)=ktop1
                           qdif(j)=qa(j,kx)-qthr0
                        end if
                    end if
                end if
            end do
        end subroutine diagnose_convection
end module phy_convmf
