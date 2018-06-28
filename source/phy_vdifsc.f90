module phy_vdifsc

    use rp_emulator
    use mod_prec

    implicit none

    private
    public vdifsc, setup_vertical_diffusion, truncate_vdifsc

    namelist /vertical_diffusion/ trshc, trvdi, trvds, redshc, rhgrad, segrad

    ! Relaxation time (in hours) for shallow convection
    type(rpe_var) :: trshc

    ! Relaxation time (in hours) for moisture diffusion
    type(rpe_var) :: trvdi

    ! Relaxation time (in hours) for super-adiab. conditions
    type(rpe_var) :: trvds

    ! Reduction factor of shallow conv. in areas of deep conv.
    type(rpe_var) :: redshc

    ! Maximum gradient of relative humidity (d_RH/d_sigma)
    type(rpe_var) :: rhgrad

    ! Minimum gradient of dry static energy (d_DSE/d_phi)
    type(rpe_var) :: segrad

    contains
        subroutine setup_vertical_diffusion(fid)
            integer, intent(in) :: fid

            read(fid, vertical_diffusion)

            write(*, vertical_diffusion)
        end subroutine setup_vertical_diffusion
        
        subroutine truncate_vdifsc()
            call apply_truncation(trshc)        
            call apply_truncation(trvdi)        
            call apply_truncation(trvds)        
            call apply_truncation(redshc)        
            call apply_truncation(rhgrad)        
            call apply_truncation(segrad)            
        end subroutine truncate_vdifsc

        subroutine vdifsc(ua,va,se,rh,qa,qsat,phi,icnv, &
                utenvd,vtenvd,ttenvd,qtenvd)
            !   subroutine vdifsc (ua,va,se,rh,qa,qsat,phi,icnv,
            !  &                   utenvd,vtenvd,ttenvd,qtenvd)
            !
            !   Purpose: Compute tendencies of momentum, energy and moisture
            !            due to vertical diffusion and shallow convection
            !   Input:   ua     = u-wind                           (3-dim)
            !            va     = v-wind                           (3-dim)
            !            se     = dry static energy                (3-dim)
            !            rh     = relative humidity [0-1]          (3-dim)
            !            qa     = specific humidity [g/kg]         (3-dim)
            !            qsat   = saturation sp. humidity [g/kg]   (3-dim)
            !            phi    = geopotential                     (3-dim)
            !            icnv   = index of deep convection         (2-dim)
            !   Output:  utenvd = u-wind tendency                  (3-dim)
            !            vtenvd = v-wind tendency                  (3-dim)
            !            ttenvd = temperature tendency             (3-dim)
            !            qtenvd = sp. humidity tendency [g/(kg s)] (3-dim)
            !

            use mod_atparam
            use mod_physcon, only: cp, alhc, sig, sigh, dsig
        
            type(rpe_var), dimension(ngp,kx), intent(in) :: &
                    ua, va, se, rh, qa, qsat, phi
            integer, intent(in) :: icnv(ngp)
            type(rpe_var), dimension(ngp,kx), intent(inout) :: &
                    utenvd, vtenvd, ttenvd, qtenvd

            integer :: nl1, j, k, k1
            type(rpe_var) :: cshc, cvdi, fshcq, fshcse, fvdiq, fvdise, drh0, &
                    fvdiq2, dmse, drh, fluxse, fluxq, fcnv, se0, one
            type(rpe_var), dimension(kx) :: rsig, rsig1

            one = 1.0_dp

            ! 1. Initalization

            ! N.B. In this routine, fluxes of dry static energy and humidity
            !      are scaled in such a way that:
            !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma

            nl1  = kx-1
            cshc = dsig(kx)/rpe_literal(3600.0_dp)
            cvdi = (sigh(nl1)-sigh(1))/((nl1-1)*rpe_literal(3600.0_dp))

            fshcq  = cshc/trshc
            fshcse = cshc/(trshc*cp)

            fvdiq  = cvdi/trvdi
            fvdise = cvdi/(trvds*cp)

            do k=1,nl1
                rsig(k)=one/dsig(k)
                rsig1(k)=one/(one-sigh(k))
            end do
            rsig(kx)=one/dsig(kx)

            utenvd = 0.0_dp
            vtenvd = 0.0_dp
            ttenvd = 0.0_dp
            qtenvd = 0.0_dp

            ! 2. Shallow convection
            drh0   = rhgrad*(sig(kx)-sig(nl1))
            fvdiq2 = fvdiq*sigh(nl1)

            do j=1,ngp
                dmse = (se(j,kx)-se(j,nl1))+alhc*(qa(j,kx)-qsat(j,nl1))
                drh  = rh(j,kx)-rh(j,nl1)
                fcnv = 1.0_dp

                if (dmse.ge.rpe_literal(0.0_dp)) then
                    if (icnv(j).gt.0) fcnv = redshc

                    fluxse         = fcnv*fshcse*dmse
                    ttenvd(j,nl1)  = fluxse*rsig(nl1)
                    ttenvd(j,kx) =-fluxse*rsig(kx)

                    if (drh.ge.rpe_literal(0.0_dp)) then
                        fluxq          = fcnv*fshcq*qsat(j,kx)*drh
                        qtenvd(j,nl1)  = fluxq*rsig(nl1)
                        qtenvd(j,kx) =-fluxq*rsig(kx)
                    end if
                else if (drh.ge.drh0) then
                  fluxq          = fvdiq2*qsat(j,nl1)*drh
                  qtenvd(j,nl1)  = fluxq*rsig(nl1)
                  qtenvd(j,kx) =-fluxq*rsig(kx)
                end if
            end do

            ! 3. Vertical diffusion of moisture above the PBL
            do k=3,kx-2
                if (sigh(k).gt.rpe_literal(0.5_dp)) then
                    drh0   = rhgrad*(sig(k+1)-sig(k))
                    fvdiq2 = fvdiq*sigh(k)

                    do j=1,ngp
                        drh=rh(j,k+1)-rh(j,k)
                        if (drh.ge.drh0) then
                            fluxq        = fvdiq2*qsat(j,k)*drh
                            qtenvd(j,k)  = qtenvd(j,k)  +fluxq*rsig(k)
                            qtenvd(j,k+1)= qtenvd(j,k+1)-fluxq*rsig(k+1)
                        end if
                    end do
                end if
            end do

            ! 4. Damping of super-adiabatic lapse rate
            do k=1,nl1
                do j=1,ngp
                    se0 = se(j,k+1)+segrad*(phi(j,k)-phi(j,k+1))

                    if (se(j,k).lt.se0) then
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
