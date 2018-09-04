module phy_vdifsc
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public vdifsc, setup_vertical_diffusion, ini_vdifsc

    ! Variables loaded in by namelist
    namelist /vertical_diffusion/ trshc, trvdi, trvds, redshc, rhgrad, segrad

    ! Relaxation time (in hours) for shallow convection
    real(dp) :: trshc

    ! Relaxation time (in hours) for moisture diffusion
    real(dp) :: trvdi

    ! Relaxation time (in hours) for super-adiab. conditions
    real(dp) :: trvds

    ! Reduction factor of shallow conv. in areas of deep conv.
    real(dp) :: redshc

    ! Maximum gradient of relative humidity (d_RH/d_sigma)
    real(dp) :: rhgrad

    ! Minimum gradient of dry static energy (d_DSE/d_phi)
    real(dp) :: segrad

    ! Local derived variables
    real(dp) :: fshcq, fshcse, fvdiq, fvdise
    real(dp), allocatable :: rsig(:), rsig1(:), drh0(:), fvdiq2(:)

    contains
        subroutine setup_vertical_diffusion(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, vertical_diffusion)

            write(*, vertical_diffusion)

            allocate(rsig(kx))
            allocate(rsig1(kx))
            allocate(drh0(3:kxm))
            allocate(fvdiq2(3:kxm))
        end subroutine setup_vertical_diffusion

        subroutine ini_vdifsc()
            ! Calculate local variables for turbulence scheme
            use mod_physcon, only: cp, alhc, sig, sigh, dsig

            real(dp) :: cshc, cvdi

            cshc = dsig(kx)/3600.0_dp
            cvdi = (sigh(kxm)-sigh(1))/((kxm-1)*3600.0_dp)

            fshcq  = cshc/trshc
            fshcse = cshc/(trshc*cp)
            fvdise = cvdi/(trvds*cp)

            rsig=1.0_dp/dsig
            rsig1=1.0_dp/(1.0_dp-sigh)

            drh0(3:kxm)   = rhgrad*(sig(4:kx)-sig(3:kxm))
            fvdiq2(3:kxm) = (cvdi/trvdi)*sigh(3:kxm)

        end subroutine ini_vdifsc

        subroutine vdifsc(ua, va, se, rh, qa, qsat, phi, icnv, &
                utenvd, vtenvd, ttenvd, qtenvd)
            !   subroutine vdifsc (ua,va,se,rh,qa,qsat,phi,icnv,
            !  &                   utenvd,vtenvd,ttenvd,qtenvd)
            !
            !   Purpose: Compute tendencies of momentum, energy and moisture
            !            due to vertical diffusion and shallow convection
            use mod_physcon, only: alhc, sigh

            !   Input:   ua     = u-wind                           (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: ua
            !            va     = v-wind                           (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: va
            !            se     = dry static energy                (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: se
            !            rh     = relative humidity [0-1]          (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: rh
            !            qa     = specific humidity [g/kg]         (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: qa
            !            qsat   = saturation sp. humidity [g/kg]   (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: qsat
            !            phi    = geopotential                     (3-dim)
            real(dp), dimension(ngp,kx), intent(in) :: phi
            !            icnv   = index of deep convection         (2-dim)
            integer, intent(in) :: icnv(ngp)

            !   Output:  utenvd = u-wind tendency                  (3-dim)
            real(dp), dimension(ngp,kx), intent(out) :: utenvd
            !            vtenvd = v-wind tendency                  (3-dim)
            real(dp), dimension(ngp,kx), intent(out) :: vtenvd
            !            ttenvd = temperature tendency             (3-dim)
            real(dp), dimension(ngp,kx), intent(out) :: ttenvd
            !            qtenvd = sp. humidity tendency [g/(kg s)] (3-dim)
            real(dp), dimension(ngp,kx), intent(out) :: qtenvd

            ! Local variables
            integer :: j, k, k1
            real(dp) :: dmse, drh, fluxse, fluxq, fcnv, se0

            ! 1. Initalization
            ! N.B. In this routine, fluxes of dry static energy and humidity
            !      are scaled in such a way that:
            !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma
            utenvd = 0.0_dp
            vtenvd = 0.0_dp
            ttenvd = 0.0_dp
            qtenvd = 0.0_dp

            ! 2. Shallow convection
            do j=1,ngp
                dmse = (se(j,kx)-se(j,kxm))+alhc*(qa(j,kx)-qsat(j,kxm))
                drh  = rh(j,kx)-rh(j,kxm)
                fcnv = 1.0_dp

                if (dmse >= 0.0_dp) then
                    if (icnv(j) > 0) fcnv = redshc

                    fluxse         = fcnv*fshcse*dmse
                    ttenvd(j,kxm)  = fluxse*rsig(kxm)
                    ttenvd(j,kx) =-fluxse*rsig(kx)

                    if (drh >= 0.0_dp) then
                        fluxq          = fcnv*fshcq*qsat(j,kx)*drh
                        qtenvd(j,kxm)  = fluxq*rsig(kxm)
                        qtenvd(j,kx) =-fluxq*rsig(kx)
                    end if
                else if (drh >= drh0(kxm)) then
                    fluxq         = fvdiq2(kxm)*qsat(j,kxm)*drh
                    qtenvd(j,kxm) = fluxq*rsig(kxm)
                    qtenvd(j,kx)  =-fluxq*rsig(kx)
                end if
            end do

            ! 3. Vertical diffusion of moisture above the PBL
            do k=3,kx-2
                if (sigh(k) > 0.5_dp) then
                    do j=1,ngp
                        drh=rh(j,k+1)-rh(j,k)
                        if (drh >= drh0(k)) then
                            fluxq        = fvdiq2(k)*qsat(j,k)*drh
                            qtenvd(j,k)  = qtenvd(j,k)  +fluxq*rsig(k)
                            qtenvd(j,k+1)= qtenvd(j,k+1)-fluxq*rsig(k+1)
                        end if
                    end do
                end if
            end do

            ! 4. Damping of super-adiabatic lapse rate
            do k=1,kxm
                do j=1,ngp
                    se0 = se(j,k+1)+segrad*(phi(j,k)-phi(j,k+1))

                    if (se(j,k) < se0) then
                        fluxse      = fvdise*(se0-se(j,k))
                        ttenvd(j,k) = ttenvd(j,k)+fluxse*rsig(k)
                        do k1=k+1,kx
                            ttenvd(j,k1) = ttenvd(j,k1)-fluxse*rsig1(k)
                        end do
                    end if
                end do
            end do
        end subroutine vdifsc
end module phy_vdifsc
