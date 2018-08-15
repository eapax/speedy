module phy_convmf
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public convmf, setup_convection, ini_convmf, truncate_convmf

    ! Variables loaded in by namelist
    namelist /convection/ psmin, trcnv, rhbl, rhil, entmax, smf

    ! Minimum (norm.) sfc. pressure for the occurrence of convection
    type(rpe_var) :: psmin

    ! Time of relaxation (in hours) towards reference state
    real(dp) :: trcnv

    ! Relative hum. threshold in the boundary layer
    type(rpe_var) :: rhbl

    ! Rel. hum. threshold in intermed. layers for secondary mass flux
    type(rpe_var) :: rhil

    ! Max. entrainment as a fraction of cloud-base mass flux
    type(rpe_var) :: entmax

    ! Ratio between secondary and primary mass flux at cloud-base
    type(rpe_var) :: smf

    ! Local derived variables
    type(rpe_var) :: fm0

    ! Local copies of mod_atparam variables
    type(rpe_var) :: alhc_cnv
    type(rpe_var), allocatable :: sig_cnv(:), wvi_cnv(:,:)

    contains
        subroutine setup_convection(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, convection)

            write(*, convection)

            allocate(sig_cnv(kx))
            allocate(wvi_cnv(kx,2))
        end subroutine setup_convection

        subroutine ini_convmf()
            ! Calculate local variables for convection scheme
            use mod_physcon, only: p0, gg, alhc, sig, dsig, wvi

            fm0=p0*dsig(kx)/(gg*trcnv*3600.0_dp)

            alhc_cnv = alhc
            sig_cnv = sig
            wvi_cnv = wvi
        end subroutine ini_convmf

        subroutine truncate_convmf()
            ! Truncate local variables for convection scheme
            ! Namelist variables
            call apply_truncation(psmin)
            call apply_truncation(rhbl)
            call apply_truncation(rhil)
            call apply_truncation(entmax)
            call apply_truncation(smf)

            ! Derived variables
            call apply_truncation(fm0)

            ! Local copies of mod_physcon
            call apply_truncation(alhc_cnv)
            call apply_truncation(sig_cnv)
            call apply_truncation(wvi_cnv)
        end subroutine truncate_convmf

        subroutine convmf (&
                psa_in,se_in,qa_in,qsat_in, hflx2tend_in, flx2tend_in, &
                itop,cbmf,precnv,dfse,dfqa)
            ! SUBROUTINE CONVMF (PSA,SE,QA,QSAT, ITOP,CBMF,PRECNV,DFSE,DFQA)
            !
            ! Purpose: Compute convective fluxes of dry static energy and
            !          moisture using a simplified mass-flux scheme

            ! Input:  PSA    = norm. surface pressure [p/p0]            (2-dim)
            type(rpe_var), intent(in) :: psa_in(ngp)
            !         SE     = dry static energy                        (3-dim)
            type(rpe_var), intent(in) :: se_in(ngp, kx)
            !         QA     = specific humidity [g/kg]                 (3-dim)
            type(rpe_var), intent(in) :: qa_in(ngp, kx)
            !         QSAT   = saturation spec. hum. [g/kg]             (3-dim)
            type(rpe_var), intent(in) :: qsat_in(ngp, kx)
            !         hflx2tend = Conversion factor between heat fluxes and T tendency
            type(rpe_var), intent(in) :: hflx2tend_in(ngp,kx)
            !         flx2tend = Conversion factor between fluxes and tendencies
            type(rpe_var), intent(in) :: flx2tend_in(ngp,kx)

            ! Output: ITOP   = top of convection (layer index)          (2-dim)
            integer, intent(out) :: itop(ngp)
            !         CBMF   = cloud-base mass flux                     (2-dim)
            type(rpe_var), intent(out) :: cbmf(ngp)
            !         PRECNV = convective precipitation [g/(m^2 s)]     (2-dim)
            type(rpe_var), intent(out) :: precnv(ngp)
            !         DFSE   = net flux of d.s.en. into each atm. layer (3-dim)
            type(rpe_var), intent(out) :: dfse(ngp,kx)
            !         DFQA   = net flux of sp.hum. into each atm. layer (3-dim)
            type(rpe_var), intent(out) :: dfqa(ngp,kx)

            ! Local copies of input variables (to be truncated)
            type(rpe_var) :: psa(ngp), se(ngp,kx), qa(ngp,kx), qsat(ngp,kx), &
                    hflx2tend(ngp,kx), flx2tend(ngp,kx)

            ! Local variables
            integer :: j, k, k1, ktop1, ktop2, nl1, nlp
            type(rpe_var) :: mss(ngp,2:kx), mse0, mse1, mss0, mss2, msthr, &
                    qdif(ngp), entr(2:kx-1), delq, enmass, fdq, fds, &
                    fmass, fpsa, fqmax, fsq, fuq, fus, qb, qmax, qsatb, qthr0, &
                    qthr1, rdps, rlhc, sb, sentr
            logical :: lqthr

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            psa = psa_in
            se = se_in
            qa = qa_in
            qsat = qsat_in
            hflx2tend = hflx2tend_in
            flx2tend = flx2tend_in

            ! 1. Initialization of output and workspace arrays
            nl1=kx-1
            nlp=kx+1
            fqmax=5.0_dp

            rdps=rpe_literal(2.0_dp)/(rpe_literal(1.0_dp)-psmin)

            dfse = 0.0_dp
            dfqa = 0.0_dp

            cbmf = 0.0_dp
            precnv = 0.0_dp

            ! Saturation moist static energy
            do k=2,kx
                do j=1,ngp
                    mss(j,k)=se(j,k)+alhc_cnv*qsat(j,k)
                end do
            end do

            ! Entrainment profile (up to sigma = 0.5)
            sentr=0.0_dp
            do k=2,nl1
                entr(k)=(max(rpe_literal(0.0_dp), &
                        sig_cnv(k)-rpe_literal(0.5_dp)))**2
                sentr=sentr+entr(k)
            end do

            sentr=entmax/sentr
            entr(2:nl1) = entr(2:nl1) * sentr

            ! 2. Check of conditions for convection
            rlhc=rpe_literal(1.0_dp)/alhc_cnv

            do j=1,ngp
                itop(j)=nlp

                if (psa(j).gt.psmin) then
                    ! Minimum of moist static energy in the lowest two levels
                    mse0=se(j,kx)+alhc_cnv*qa(j,kx)
                    mse1=se(j,nl1) +alhc_cnv*qa(j,nl1)
                    mse1=min(mse0,mse1)

                    ! Saturation (or super-saturated) moist static energy in PBL
                    mss0=max(mse0,mss(j,kx))

                    ktop1=kx
                    ktop2=kx

                    do k=kx-3,3,-1
                        mss2=mss(j,k)+wvi_cnv(k,2)*(mss(j,k+1)-mss(j,k))

                        ! Check 1: conditional instability
                        !          (MSS in PBL > MSS at top level)
                        if (mss0.gt.mss2) then
                           ktop1=k
                        end if

                        ! Check 2: gradient of actual moist static energy
                        !          between lower and upper troposphere
                        if (mse1.gt.mss2) then
                           ktop2=k
                           msthr=mss2
                        end if
                    end do

                    if (ktop1.lt.kx) then
                        ! Check 3: RH > RH_c at both k=kx and k=NL1
                        qthr0=rhbl*qsat(j,kx)
                        qthr1=rhbl*qsat(j,nl1)
                        lqthr=(qa(j,kx).gt.qthr0.and.qa(j,nl1).gt.qthr1)

                        if (ktop2.lt.kx) then
                           itop(j)=ktop1
                           qdif(j)=max(qa(j,kx)-qthr0,(mse0-msthr)*rlhc)
                        else if (lqthr) then
                           itop(j)=ktop1
                           qdif(j)=qa(j,kx)-qthr0
                        end if
                    end if
                end if
            end do

            ! 3. Convection over selected grid-points
            do j=1,ngp
                if (itop(j).eq.nlp) cycle

                ! 3.1 Boundary layer (cloud base)
                k =kx
                k1=k-1

                ! Maximum specific humidity in the PBL
                qmax=max(rpe_literal(1.01_dp)*qa(j,k),qsat(j,k))

                ! Dry static energy and moisture at upper boundary
                sb=se(j,k1)+wvi_cnv(k1,2)*(se(j,k)-se(j,k1))
                qb=qa(j,k1)+wvi_cnv(k1,2)*(qa(j,k)-qa(j,k1))
                qb=min(qb,qa(j,k))

                ! Cloud-base mass flux, computed to satisfy:
                ! fmass*(qmax-qb)*(g/dp)=qdif/trcnv
                fpsa=psa(j)*min(rpe_literal(1.0_dp),(psa(j)-psmin)*rdps)
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
                    sb=se(j,k1)+wvi_cnv(k1,2)*(se(j,k)-se(j,k1))
                    qb=qa(j,k1)+wvi_cnv(k1,2)*(qa(j,k)-qa(j,k1))
                    fds=fmass*sb
                    fdq=fmass*qb

                    ! Net flux of dry static energy and moisture
                    dfse(j,k)=dfse(j,k)+fds-fus
                    dfqa(j,k)=dfqa(j,k)+fdq-fuq

                    ! Secondary moisture flux
                    delq=rhil*qsat(j,k)-qa(j,k)
                    if (delq.gt.rpe_literal(0.0_dp)) then
                        fsq=smf*cbmf(j)*delq
                        dfqa(j,k)   =dfqa(j,k)   +fsq
                        dfqa(j,kx)=dfqa(j,kx)-fsq
                    end if
                end do

                ! 3.3 Top layer (condensation and detrainment)
                k=itop(j)

                ! Flux of convective precipitation
                qsatb=qsat(j,k)+wvi_cnv(k,2)*(qsat(j,k+1)-qsat(j,k))
                precnv(j)=max(fuq-fmass*qsatb,rpe_literal(0.0_dp))

                ! Net flux of dry static energy and moisture
                dfse(j,k)=fus-fds+alhc_cnv*precnv(j)
                dfqa(j,k)=fuq-fdq-precnv(j)
            end do

            ! Convert fluxes to temperature tendencies
            dfse = dfse*hflx2tend
            dfqa = dfqa*flx2tend
        end subroutine convmf
end module phy_convmf
