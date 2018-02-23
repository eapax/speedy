module vertical_diffusion
    implicit none

    private
    public vdifsc

    namelist /vertical_diffusion/ trshc, trvdi, trvds, redshc, rhgrad, segrad

    ! Relaxation time (in hours) for shallow convection
    real, parameter :: trshc = 6.0

    ! Relaxation time (in hours) for moisture diffusion
    real, parameter :: trvdi = 24.0

    ! Relaxation time (in hours) for super-adiab. conditions
    real, parameter :: trvds = 6.0

    ! Reduction factor of shallow conv. in areas of deep conv.
    real, parameter :: redshc = 0.5

    ! Maximum gradient of relative humidity (d_RH/d_sigma)
    real, parameter :: rhgrad = 0.5

    ! Minimum gradient of dry static energy (d_DSE/d_phi)
    real, parameter :: segrad = 0.1

    contains
        subroutine setup_vertical_diffusion(fid)
            integer, intent(in) :: fid

            read(fid, vertical_diffusion)
        end subroutine setup_vertical_diffusion

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
        
            integer, parameter :: ngp=ix*il
        
            real, dimension(ngp,kx), intent(in) :: ua, va, se, rh, qa, qsat, phi
            integer, intent(in) :: icnv(ngp)
            real, dimension(ngp,kx), intent(inout) :: &
                    utenvd, vtenvd, ttenvd, qtenvd
        
            integer :: nl1, j, k, k1
            real :: cshc, cvdi, fshcq, fshcse, fvdiq, fvdise, drh0, fvdiq2, &
                    dmse, drh, fluxse, fluxq, fcnv, se0
            real, dimension(kx) :: rsig, rsig1
        
            ! 1. Initalization
        
            ! N.B. In this routine, fluxes of dry static energy and humidity
            !      are scaled in such a way that:
            !      d_T/dt = d_F'(SE)/d_sigma,  d_Q/dt = d_F'(Q)/d_sigma
        
            nl1  = kx-1
            cshc = dsig(kx)/3600.
            cvdi = (sigh(nl1)-sigh(1))/((nl1-1)*3600.)
        
            fshcq  = cshc/trshc
            fshcse = cshc/(trshc*cp)
        
            fvdiq  = cvdi/trvdi
            fvdise = cvdi/(trvds*cp)
        
            do k=1,nl1
                rsig(k)=1./dsig(k)
                rsig1(k)=1./(1.-sigh(k))
            end do
            rsig(kx)=1./dsig(kx)
        
            utenvd = 0.0
            vtenvd = 0.0
            ttenvd = 0.0
            qtenvd = 0.0
           
            ! 2. Shallow convection
            drh0   = rhgrad*(sig(kx)-sig(nl1))
            fvdiq2 = fvdiq*sigh(nl1)
        
            do j=1,ngp
                dmse = (se(j,kx)-se(j,nl1))+alhc*(qa(j,kx)-qsat(j,nl1))
                drh  = rh(j,kx)-rh(j,nl1)
                fcnv = 1.
        
                if (dmse.ge.0.0) then
                    if (icnv(j).gt.0) fcnv = redshc
          
                    fluxse         = fcnv*fshcse*dmse
                    ttenvd(j,nl1)  = fluxse*rsig(nl1)
                    ttenvd(j,kx) =-fluxse*rsig(kx)
          
                    if (drh.ge.0.0) then
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
                if (sigh(k).gt.0.5) then
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
end module vertical_diffusion
