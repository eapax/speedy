module phy_lscond

    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public lscond, setup_condensation, truncate_lscond

    namelist /condensation/ trlsc, rhlsc, drhlsc, rhblsc

    ! Relaxation time (in hours) for specific humidity 
    type(rpe_var) :: trlsc

    ! Maximum relative humidity threshold (at sigma=1)
    type(rpe_var) :: rhlsc

    ! Vertical range of relative humidity threshold
    type(rpe_var) :: drhlsc

    ! Relative humidity threshold for boundary layer
    type(rpe_var) :: rhblsc
    
    contains
        subroutine setup_condensation(fid)
            integer, intent(in) :: fid

            read(fid, condensation)

            write(*, condensation)
        end subroutine setup_condensation
        
        subroutine truncate_lscond()
            call apply_truncation(trlsc)        
            call apply_truncation(rhlsc)        
            call apply_truncation(drhlsc)        
            call apply_truncation(rhblsc)            
        end subroutine truncate_lscond

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

            type(rpe_var), intent(in) :: psa(ngp), qa(ngp,kx), qsat(ngp,kx)
        
            integer, intent(inout) :: itop(ngp)
            type(rpe_var), intent(inout) :: precls(ngp), dtlsc(ngp,kx), dqlsc(ngp,kx)

            integer :: j, k
            type(rpe_var) :: psa2(ngp), dqa, dqmax, pfact, prg, qsmax, rhref, rtlsc, sig2, tfact

            ! 1. Initialization
            qsmax = 10.0_dp

            rtlsc = rpe_literal(1.0_dp)/(trlsc*rpe_literal(3600.0_dp))
            tfact = alhc/cp
            prg = p0/gg

            dtlsc(:,1) = 0.0_dp
            dqlsc(:,1) = 0.0_dp
            precls  = 0.0_dp
            do j=1,ngp
                psa2(j) = psa(j)*psa(j)
            end do

            ! 2. Tendencies of temperature and moisture
            !    NB. A maximum heating rate is imposed to avoid
            !        grid-point-storm instability
            do k=2,kx
                sig2=sig(k)*sig(k)
                rhref = rhlsc+drhlsc*(sig2-rpe_literal(1.0_dp))
                if (k.eq.kx) rhref = max(rhref,rhblsc)
                dqmax = qsmax*sig2*rtlsc

                do j=1,ngp
                    dqa = rhref*qsat(j,k)-qa(j,k)
                    if (dqa.lt.rpe_literal(0.0_dp)) then
                        itop(j)    = min(k,itop(j))
                        dqlsc(j,k) = dqa*rtlsc
                        dtlsc(j,k) = tfact*min(-dqlsc(j,k),dqmax*psa2(j))
                    else
                        dqlsc(j,k) = 0.0_dp
                        dtlsc(j,k) = 0.0_dp
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
