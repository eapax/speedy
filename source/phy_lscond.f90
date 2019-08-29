module phy_lscond
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public lscond, setup_condensation, ini_lscond, truncate_lscond

    ! Variables loaded in by namelist
    namelist /condensation/ trlsc, rhlsc, drhlsc, rhblsc

    real(dp), parameter :: qsmax = 10.0_dp

    ! Relaxation time (in hours) for specific humidity
    real(dp) :: trlsc

    ! Maximum relative humidity threshold (at sigma=1)
    real(dp) :: rhlsc

    ! Vertical range of relative humidity threshold
    real(dp) :: drhlsc

    ! Relative humidity threshold for boundary layer
    real(dp) :: rhblsc

    ! Local derived variables
    type(rpe_var) :: rtlsc, tfact
    type(rpe_var), allocatable :: pfact(:), rhref(:), dqmax(:)

    contains
        subroutine setup_condensation(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, condensation)

            write(*, condensation)

            allocate(pfact(kx))
            allocate(rhref(2:kx))
            allocate(dqmax(2:kx))
        end subroutine setup_condensation

        subroutine ini_lscond()
            ! Calculate local variables for condensation scheme
            use mod_physcon, only: p0, gg, cp, alhc, sig, dsig

            real(dp) :: sig2
            integer :: k

            rtlsc = 1.0_dp/trlsc
            tfact = alhc/cp
            pfact = dsig*p0/gg

            do k=2,kx
                sig2=sig(k)*sig(k)
                rhref(k) = rhlsc+drhlsc*(sig2-1.0_dp)
                dqmax(k) = qsmax*sig2*rtlsc
            end do
            rhref(kx) = max(rhref(kx)%val, rhblsc)

        end subroutine ini_lscond

        subroutine truncate_lscond()
            ! Truncate local variables for condensation scheme
            ! Derived variables
            call apply_truncation(rtlsc)
            call apply_truncation(tfact)
            call apply_truncation(pfact)
            call apply_truncation(rhref)
            call apply_truncation(dqmax)
        end subroutine truncate_lscond

        subroutine lscond(psa_in, qa_in, qsat_in, itop, &
                precls_out, dtlsc_out, dqlsc_out)
            !  subroutine lscond (psa,qa,qsat,itop,precls,dtlsc,dqlsc)
            !
            !  Purpose: Compute large-scale precipitation and
            !           associated tendencies of temperature and moisture
            !  Input:   psa    = norm. surface pressure [p/p0]           (2-dim)
            real(dp), intent(in) :: psa_in(ngp)
            !           qa     = specific humidity [g/kg]                (3-dim)
            real(dp), intent(in) :: qa_in(ngp,kx)
            !           qsat   = saturation spec. hum. [g/kg]            (3-dim)
            real(dp), intent(in) :: qsat_in(ngp,kx)

            !           itop   = top of convection (layer index)         (2-dim)
            !  Output:  itop   = top of conv+l.s.condensat.(layer index) (2-dim)
            integer, intent(inout) :: itop(ngp)
            !           precls = large-scale precipitation [g/(m^2 s)]   (2-dim)
            real(dp), intent(out) :: precls_out(ngp)
            !           dtlsc  = temperature tendency from l.s. cond     (3-dim)
            real(dp), intent(out) :: dtlsc_out(ngp,kx)
            !           dqlsc  = hum. tendency [g/(kg s)] from l.s. cond (3-dim)
            real(dp), intent(out) :: dqlsc_out(ngp,kx)

            ! Local copies of input variables
            type(rpe_var) :: psa(ngp), qa(ngp,kx), qsat(ngp,kx)

            ! Local copies of output variables
            type(rpe_var) :: precls(ngp), dtlsc(ngp,kx), dqlsc(ngp,kx)

            ! Local variables
            integer :: j, k
            type(rpe_var) ::  dqa

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            psa = psa_in
            qa = qa_in
            qsat = qsat_in

            ! 1. Initialization
            dtlsc(:,1) = 0.0_dp
            dqlsc(:,1) = 0.0_dp
            precls  = 0.0_dp

            ! 2. Tendencies of temperature and moisture
            !    NB. A maximum heating rate is imposed to avoid
            !        grid-point-storm instability
            do k=2,kx
                do j=1,ngp
                    dqa = rhref(k)*qsat(j,k)-qa(j,k)
                    if (dqa<rpe_literal(0.0_dp)) then
                        itop(j)    = min(k,itop(j))
                        dqlsc(j,k) = dqa*rtlsc
                        dtlsc(j,k) = tfact*min(-dqlsc(j,k), dqmax(k)*psa(j)**2)
                    else
                        dqlsc(j,k) = 0.0_dp
                        dtlsc(j,k) = 0.0_dp
                    endif
                end do
            end do

            ! 3. Large-scale precipitation
            do k=2,kx
                do j=1,ngp
                    precls(j) = precls(j)-pfact(k)*dqlsc(j,k)
                end do
            end do

            do j=1,ngp
                precls(j) = precls(j)*psa(j)
            end do

            precls_out = precls
            dtlsc_out = dtlsc
            dqlsc_out = dqlsc
        end subroutine lscond
end module phy_lscond
