module phy_lscond
    implicit none

    private
    public lscond, setup_condensation

    namelist /condensation/ trlsc, rhlsc, drhlsc, rhblsc

    ! Relaxation time (in hours) for specific humidity 
    real :: trlsc  = 4.0

    ! Maximum relative humidity threshold (at sigma=1)
    real :: rhlsc  = 0.9

    ! Vertical range of relative humidity threshold
    real :: drhlsc = 0.1

    ! Relative humidity threshold for boundary layer
    real :: rhblsc = 0.95
    
    contains
        subroutine setup_condensation(fid)
            integer, intent(in) :: fid

            read(fid, condensation)

            write(*, condensation)
        end subroutine setup_condensation
    
        subroutine lscond(psa,qa,qsat,itop,precls,dtlsc,dqlsc)
            !  subroutine lscond (psa,qa,qsat,itop,precls,dtlsc,dqlsc)
            !
            !  Purpose: Compute large-scale precipitation and
            !           associated tendencies of temperature and moisture
            !  Input:   psa    = norm. surface pressure [p/p0]           (2-dim)
            !           qa     = specific humidity [g/kg]                (3-dim)
            !           qsat   = saturation spec. hum. [g/kg]            (3-dim)
            !           itop   = top of convection (layer index)         (2-dim)
            !  Output:  itop   = top of conv+l.s.condensat.(layer index) (2-dim)
            !           precls = large-scale precipitation [g/(m^2 s)]   (2-dim)
            !           dtlsc  = temperature tendency from l.s. cond     (3-dim)
            !           dqlsc  = hum. tendency [g/(kg s)] from l.s. cond (3-dim)
        
            use mod_atparam
            use mod_physcon, only: p0, gg, cp, alhc, alhs, sig, dsig
        
            real, intent(in) :: psa(ngp), qa(ngp,kx), qsat(ngp,kx)
        
            integer, intent(inout) :: itop(ngp)
            real, intent(inout) :: precls(ngp), dtlsc(ngp,kx), dqlsc(ngp,kx)
        
            integer :: j, k
            real :: psa2(ngp), dqa, dqmax, pfact, prg, qsmax, rhref, rtlsc, sig2, tfact
        
            ! 1. Initialization
            qsmax = 10.
        
            rtlsc = 1./(trlsc*3600.)
            tfact = alhc/cp
            prg = p0/gg
        
            dtlsc(:,1) = 0.
            dqlsc(:,1) = 0.
            precls  = 0.
            do j=1,ngp
                psa2(j) = psa(j)*psa(j)
            end do
        
            ! 2. Tendencies of temperature and moisture
            !    NB. A maximum heating rate is imposed to avoid 
            !        grid-point-storm instability 
            do k=2,kx
                sig2=sig(k)*sig(k)
                rhref = rhlsc+drhlsc*(sig2-1.)
                if (k.eq.kx) rhref = max(rhref,rhblsc)
                dqmax = qsmax*sig2*rtlsc
        
                do j=1,ngp
                    dqa = rhref*qsat(j,k)-qa(j,k)
                    if (dqa.lt.0.0) then
                        itop(j)    = min(k,itop(j))
                        dqlsc(j,k) = dqa*rtlsc
                        dtlsc(j,k) = tfact*min(-dqlsc(j,k),dqmax*psa2(j))
                    else
                        dqlsc(j,k) = 0.
                        dtlsc(j,k) = 0.
                    endif
                end do
            end do
        
            ! 3. Large-scale precipitation
            do k=2,kx
                pfact = dsig(k)*prg
                do j=1,ngp
                    precls(j) = precls(j)-pfact*dqlsc(j,k)
                end do
            end do
        
            do j=1,ngp
                precls(j) = precls(j)*psa(j)
            end do
        end subroutine lscond
end module phy_lscond
