module phy_vdifsc
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    private
    public vdifsc, setup_vertical_diffusion, ini_vdifsc, truncate_vdifsc

    ! Variables loaded in by namelist
    namelist /vertical_diffusion/ trshc, trvdi, trvds, redshc, rhgrad, segrad

    ! Relaxation time (in hours) for shallow convection
    real(dp) :: trshc

    ! Relaxation time (in hours) for moisture diffusion
    real(dp) :: trvdi

    ! Relaxation time (in hours) for super-adiab. conditions
    real(dp) :: trvds

    ! Reduction factor of shallow conv. in areas of deep conv.
    type(rpe_var) :: redshc

    ! Maximum gradient of relative humidity (d_RH/d_sigma)
    type(rpe_var) :: rhgrad

    ! Minimum gradient of dry static energy (d_DSE/d_phi)
    type(rpe_var) :: segrad

    ! Local derived variables
    type(rpe_var) :: fshcq, fshcse, fvdise
    type(rpe_var), allocatable :: rsig(:), rsig1(:), drh0(:), fvdiq2(:)

    ! Local copies of mod_atparam variables
    type(rpe_var) :: alhc_vdif, cp_vdif
    type(rpe_var), allocatable :: sigh_vdif(:)

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
            allocate(sigh_vdif(0:kx))
        end subroutine setup_vertical_diffusion

        subroutine ini_vdifsc()
            ! Calculate local variables for turbulence scheme
            use mod_physcon, only: cp, alhc, sig, sigh, dsig

            real(dp) :: cshc, cvdi

            cshc = dsig(kx)
            cvdi = (sigh(kxm)-sigh(1))/(kxm-1)

            fshcq  = cshc/trshc
            fshcse = cshc/(trshc*cp)
            fvdise = cvdi/(trvds*cp)

            rsig=1.0_dp/dsig
            rsig1=1.0_dp/(1.0_dp-sigh)

            drh0(3:kxm)   = rhgrad*(sig(4:kx)-sig(3:kxm))
            fvdiq2(3:kxm) = (cvdi/trvdi)*sigh(3:kxm)

            alhc_vdif = alhc
            sigh_vdif = sigh
            cp_vdif = cp
        end subroutine ini_vdifsc

        subroutine truncate_vdifsc()
            ! Truncate local variables for turbulence scheme
            ! Namelist variables
            call apply_truncation(redshc)
            call apply_truncation(rhgrad)
            call apply_truncation(segrad)

            ! Derived variables
            call apply_truncation(fshcq)
            call apply_truncation(fshcse)
            call apply_truncation(fvdise)
            call apply_truncation(rsig)
            call apply_truncation(rsig1)
            call apply_truncation(drh0)
            call apply_truncation(fvdiq2)

            ! Local copies of mod_physcon
            call apply_truncation(alhc_vdif)
            call apply_truncation(sigh_vdif)
            call apply_truncation(cp_vdif)
        end subroutine truncate_vdifsc

        subroutine vdifsc(se_in, rh_in, qa_in, qsat_in, phi_in, icnv, &
                utenvd, vtenvd, ttenvd, qtenvd)
            !   subroutine vdifsc (se,rh,qa,qsat,phi,icnv,
            !  &                   utenvd,vtenvd,ttenvd,qtenvd)
            !
            !   Purpose: Compute tendencies of momentum, energy and moisture
            !            due to vertical diffusion and shallow convection

            !            se     = dry static energy                (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(in) :: se_in
            !            rh     = relative humidity [0-1]          (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(in) :: rh_in
            !            qa     = specific humidity [g/kg]         (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(in) :: qa_in
            !            qsat   = saturation sp. humidity [g/kg]   (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(in) :: qsat_in
            !            phi    = geopotential                     (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(in) :: phi_in
            !            icnv   = index of deep convection         (2-dim)
            integer, intent(in) :: icnv(ngp)

            !   Output:  utenvd = u-wind tendency                  (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(out) :: utenvd
            !            vtenvd = v-wind tendency                  (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(out) :: vtenvd
            !            ttenvd = temperature tendency             (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(out) :: ttenvd
            !            qtenvd = sp. humidity tendency [g/(kg s)] (3-dim)
            type(rpe_var), dimension(ngp,kx), intent(out) :: qtenvd

            ! Local copies of input variables
            type(rpe_var), dimension(ngp, kx) :: se, rh, qa, qsat, phi


            ! Local variables
            integer :: j, k, k1
            type(rpe_var) :: dmse, drh, fluxse, fluxq, fcnv, se0

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            se = se_in
            rh = rh_in
            qa = qa_in
            qsat = qsat_in
            phi = phi_in / cp_vdif

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
                dmse = (se(j,kx)-se(j,kxm))+(alhc_vdif/cp_vdif)*(qa(j,kx)-qsat(j,kxm))
                drh  = rh(j,kx)-rh(j,kxm)
                fcnv = 1.0_dp

                if (dmse>=rpe_literal(0.0_dp)) then
                    if (icnv(j)>0) fcnv = redshc

                    fluxse         = fcnv*fshcse*dmse*cp_vdif
                    ttenvd(j,kxm)  = fluxse*rsig(kxm)
                    ttenvd(j,kx) =-fluxse*rsig(kx)

                    if (drh>=rpe_literal(0.0_dp)) then
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
                if (sigh_vdif(k)>rpe_literal(0.5_dp)) then
                    do j=1,ngp
                        drh=rh(j,k+1)-rh(j,k)
                        if (drh>=drh0(k)) then
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

                    if (se(j,k)<se0) then
                        fluxse      = fvdise*(se0-se(j,k))*cp_vdif
                        ttenvd(j,k) = ttenvd(j,k)+fluxse*rsig(k)
                        do k1=k+1,kx
                            ttenvd(j,k1) = ttenvd(j,k1)-fluxse*rsig1(k)
                        end do
                    end if
                end do
            end do
        end subroutine vdifsc
end module phy_vdifsc
