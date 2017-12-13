module convection
    
    implicit none

    private
    public convmf

    ! Minimum (norm.) sfc. pressure for the occurrence of convection
    real, parameter :: psmin = 0.8

    ! Time of relaxation (in hours) towards reference state
    real, parameter :: trcnv = 6.0

    ! Relative hum. threshold in the boundary layer
    real, parameter :: rhbl = 0.9

    ! Rel. hum. threshold in intermed. layers for secondary mass flux
    real, parameter :: rhil = 0.7

    ! Max. entrainment as a fraction of cloud-base mass flux
    real, parameter :: entmax = 0.5

    ! Ratio between secondary and primary mass flux at cloud-base
    real, parameter :: smf = 0.8
    
    contains
        
        subroutine convmf (psa,se,qa,qsat,itop,cbmf,precnv,dfse,dfqa)
            ! SUBROUTINE CONVMF (PSA,SE,QA,QSAT, ITOP,CBMF,PRECNV,DFSE,DFQA)
            !
            ! Purpose: Compute convective fluxes of dry static energy and moisture
            !          using a simplified mass-flux scheme
            ! Input:  PSA    = norm. surface pressure [p/p0]            (2-dim)
            !         SE     = dry static energy                        (3-dim)
            !         QA     = specific humidity [g/kg]                 (3-dim)
            !         QSAT   = saturation spec. hum. [g/kg]             (3-dim)
            ! Output: ITOP   = top of convection (layer index)          (2-dim)
            !         CBMF   = cloud-base mass flux                     (2-dim)
            !         PRECNV = convective precipitation [g/(m^2 s)]     (2-dim)
            !         DFSE   = net flux of d.s.en. into each atm. layer (3-dim)
            !         DFQA   = net flux of sp.hum. into each atm. layer (3-dim)
            
            use mod_atparam
            use mod_physcon, only: p0, gg, alhc, alhs, sig, dsig, wvi
        
            integer, parameter :: ngp=ix*il
        
            real, intent(in) :: psa(ngp), se(ngp,kx), qa(ngp,kx), qsat(ngp,kx)
        
            integer, intent(inout) :: itop(ngp)
            real, intent(inout) :: cbmf(ngp), precnv(ngp), dfse(ngp,kx), dfqa(ngp,kx)
            !fk#if defined(KNMI)
            !fkreal :: ts(ngp),snowcv(ngp)
            !fk#endif
        
            integer :: j, k, k1, ktop1, ktop2, nl1, nlp
            real :: mss(ngp,2:kx), mse0, mse1, mss0, mss2, msthr, qdif(ngp)
            real :: entr(2:kx-1), delq, enmass, fdq, fds, fm0, fmass, fpsa, fqmax
            real :: fsq, fuq, fus, qb, qmax, qsatb, qthr0, qthr1, rdps, rlhc, sb, sentr
            logical :: lqthr
        
            ! 1. Initialization of output and workspace arrays
            nl1=kx-1
            nlp=kx+1
            fqmax=5.
        
            fm0=p0*dsig(kx)/(gg*trcnv*3600)
            rdps=2./(1.-psmin)
           
            dfse = 0.0
            dfqa = 0.0
        
            cbmf = 0.0
            precnv = 0.0
        
            ! Saturation moist static energy
            do k=2,kx
                do j=1,ngp
                    mss(j,k)=se(j,k)+alhc*qsat(j,k)
                end do
            end do
        
            ! Entrainment profile (up to sigma = 0.5)
            sentr=0.
            do k=2,nl1
                entr(k)=(max(0.,sig(k)-0.5))**2
                sentr=sentr+entr(k)
            end do
        
            sentr=entmax/sentr
            entr(2:nl1) = entr(2:nl1) * sentr
        
            ! 2. Check of conditions for convection
            rlhc=1./alhc
        
            do j=1,ngp
                itop(j)=nlp
        
                if (psa(j).gt.psmin) then
                    ! Minimum of moist static energy in the lowest two levels
                    mse0=se(j,kx)+alhc*qa(j,kx)
                    mse1=se(j,nl1) +alhc*qa(j,nl1)
                    mse1=min(mse0,mse1)
        
                    ! Saturation (or super-saturated) moist static energy in PBL 
                    mss0=max(mse0,mss(j,kx))
        
                    ktop1=kx
                    ktop2=kx
        
                    do k=kx-3,3,-1
                        mss2=mss(j,k)+wvi(k,2)*(mss(j,k+1)-mss(j,k))
            
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
                qmax=max(1.01*qa(j,k),qsat(j,k))
            
                ! Dry static energy and moisture at upper boundary
                sb=se(j,k1)+wvi(k1,2)*(se(j,k)-se(j,k1))
                qb=qa(j,k1)+wvi(k1,2)*(qa(j,k)-qa(j,k1))
                qb=min(qb,qa(j,k))
            
                ! Cloud-base mass flux, computed to satisfy:
                ! fmass*(qmax-qb)*(g/dp)=qdif/trcnv
                fpsa=psa(j)*min(1.,(psa(j)-psmin)*rdps)
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
                    if (delq.gt.0.0) then
                        fsq=smf*cbmf(j)*delq
                        dfqa(j,k)   =dfqa(j,k)   +fsq 
                        dfqa(j,kx)=dfqa(j,kx)-fsq
                    end if
                end do
        
                ! 3.3 Top layer (condensation and detrainment)
                k=itop(j)
            
                ! Flux of convective precipitation
                qsatb=qsat(j,k)+wvi(k,2)*(qsat(j,k+1)-qsat(j,k))
                precnv(j)=max(fuq-fmass*qsatb,0.0)
            
                ! Net flux of dry static energy and moisture
                dfse(j,k)=fus-fds+alhc*precnv(j)
                dfqa(j,k)=fuq-fdq-precnv(j)
            end do
        end subroutine convmf
end module convection
