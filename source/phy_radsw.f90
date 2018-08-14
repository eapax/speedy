module phy_radsw
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_sw_radiation, ini_radsw, truncate_radsw, radsw
    public nstrad, lradsw

    ! Flag for shortwave radiation routine (updated each timestep in dyn_stloop)
    logical :: lradsw

    ! Variables loaded in by namelist
    namelist /sw_radiation/ nstrad, &
            absdry, absaer, abswv1, abswv2, &
            abscl1, abscl2, &
            ablwin, ablwv1, ablwv2, &
            ablcl1, ablcl2

    ! Period (no. of steps) for shortwave radiation
    integer :: nstrad

    !          shortwave absorptivities (for dp = 10^5 Pa) :
    ! absdry = abs. of dry air      (visible band)
    type(rpe_var) :: absdry
    ! absaer = abs. of aerosols     (visible band)
    type(rpe_var) :: absaer
    ! abswv1 = abs. of water vapour (visible band, for dq = 1 g/kg)
    type(rpe_var) :: abswv1
    ! abswv2 = abs. of water vapour (near IR band, for dq = 1 g/kg)
    type(rpe_var) :: abswv2

    ! abscl2 = abs. of clouds       (visible band, for dq_base = 1 g/kg)
    type(rpe_var) :: abscl1
    ! abscl1 = abs. of clouds       (visible band, maximum value)
    type(rpe_var) :: abscl2

    !          longwave absorptivities (per dp = 10^5 Pa) :
    ! ablwin = abs. of air in "window" band
    type(rpe_var) :: ablwin
    ! ablwv1 = abs. of water vapour in H2O band 1 (weak),   for dq = 1 g/kg
    type(rpe_var) :: ablwv1
    ! ablwv2 = abs. of water vapour in H2O band 2 (strong), for dq = 1 g/kg
    type(rpe_var) :: ablwv2

    ! ablcl1 = abs. of "thick" clouds in window band (below cloud top)
    type(rpe_var) :: ablcl1
    ! ablcl2 = abs. of "thin" upper clouds in window and H2O bands
    type(rpe_var) :: ablcl2

    contains
        subroutine setup_sw_radiation(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, sw_radiation)
            write(*, sw_radiation)
        end subroutine setup_sw_radiation

        subroutine ini_radsw()
            ! Calculate local variables for short-wave radiation scheme
        end subroutine ini_radsw

        subroutine truncate_radsw()
            ! Truncate local variables for short-wave radiation scheme
            ! Namelist variables
            call apply_truncation(absdry)
            call apply_truncation(absaer)
            call apply_truncation(abswv1)
            call apply_truncation(abswv2)
            call apply_truncation(abscl1)
            call apply_truncation(abscl2)
            call apply_truncation(ablwin)
            call apply_truncation(ablwv1)
            call apply_truncation(ablwv2)
            call apply_truncation(ablcl1)
            call apply_truncation(ablcl2)
        end subroutine truncate_radsw

        subroutine radsw(psa,qa,icltop,cloudc,clstr,fsfcd,fsfc,ftop,dfabs)
            !  subroutine radsw (psa,qa,icltop,cloudc,clstr,
            ! &                  fsfcd,fsfc,ftop,dfabs)
            !
            !  purpose: compute the absorption of shortwave radiation and
            !           initialize arrays for longwave-radiation routines
            !  input:   psa    = norm. surface pressure [p/p0]           (2-dim)
            !           qa     = specific humidity [g/kg]                (3-dim)
            !           icltop = cloud top level                         (2-dim)
            !           cloudc = total cloud cover                       (2-dim)
            !           clstr  = stratiform cloud cover                  (2-dim)
            !  output:  fsfcd  = downward-only flux of sw rad. at the surface (2-dim)
            !           fsfc   = net (downw.) flux of sw rad. at the surface  (2-dim)
            !           ftop   = net (downw.) flux of sw rad. at the atm. top (2-dim)
            !           dfabs  = flux of sw rad. absorbed by each atm. layer  (3-dim)

            use mod_physcon, only: sig, dsig
            use mod_physvar, only: tau2, stratc, flux
            use mod_fordate, only: albsfc, ablco2
            use mod_solar, only: fsol, ozone, ozupp, zenit, stratz
            use phy_cloud, only: albcl, albcls, qcloud
            use phy_radlw, only: epslw

            integer, intent(in) :: icltop(ngp)
            type(rpe_var), intent(in) :: psa(ngp), qa(ngp,kx), cloudc(ngp), &
                    clstr(ngp)
            type(rpe_var), intent(inout) :: ftop(ngp), fsfc(ngp), fsfcd(ngp), &
                    dfabs(ngp,kx)

            integer :: j, k, nl1
            type(rpe_var) :: acloud(ngp), psaz(ngp), abs1, acloud1, deltap, eps1
            type(rpe_var) :: fband1, fband2

            nl1 = kx-1

            fband2 = 0.05_dp
            fband1 = rpe_literal(1.0_dp)-fband2

            ! ALBMINL=0.05
            ! ALBCLS = 0.5

            ! 1.  Initialization
            tau2 = 0.0_dp

            do j=1,ngp
                !fk-- change to ensure only icltop <= kx used
                if(icltop(j) .le. kx) then
                  tau2(j,icltop(j),3)= albcl*cloudc(j)
                endif
                !fk-- end change
                tau2(j,kx,3)     = albcls*clstr(j)
            end do

            ! 2. Shortwave transmissivity:
            ! function of layer mass, ozone (in the statosphere),
            ! abs. humidity and cloud cover (in the troposphere)

            do j=1,ngp
                psaz(j)=psa(j)*zenit(j)
                acloud(j)=cloudc(j)*min(abscl1*qcloud(j),abscl2)
            end do

            do j=1,ngp
                deltap=psaz(j)*dsig(1)
                tau2(j,1,1)=exp(-deltap*absdry)
            end do

            do k=2,nl1
                abs1=absdry+absaer*sig(k)**2
                do j=1,ngp
                    deltap=psaz(j)*dsig(k)
                    if (k.ge.icltop(j)) then
                        tau2(j,k,1)=exp(-deltap*(abs1+abswv1*qa(j,k)+acloud(j)))
                    else
                      tau2(j,k,1)=exp(-deltap*(abs1+abswv1*qa(j,k)))
                    endif
                end do
            end do

            abs1=absdry+absaer*sig(kx)**2
            do j=1,ngp
                deltap=psaz(j)*dsig(kx)
                tau2(j,kx,1)=exp(-deltap*(abs1+abswv1*qa(j,kx)))
            end do

            do k=2,kx
                do j=1,ngp
                  deltap=psaz(j)*dsig(k)
                  tau2(j,k,2)=exp(-deltap*abswv2*qa(j,k))
                end do
            end do

            ! 3. Shortwave downward flux
            ! 3.1 Initialization of fluxes
            ftop = fsol
            flux(:,1) = fsol * fband1
            flux(:,2) = fsol * fband2

            ! 3.2 Ozone and dry-air absorption in the stratosphere
            K=1
            do j=1,ngp
                dfabs(j,k)=flux(j,1)
                flux (j,1)=tau2(j,k,1)*(flux(j,1)-ozupp(j)*psa(j))
                dfabs(j,k)=dfabs(j,k)-flux(j,1)
            end do

            k=2
            do j=1,ngp
                dfabs(j,k)=flux(j,1)
                flux (j,1)=tau2(j,k,1)*(flux(j,1)-ozone(j)*psa(j))
                dfabs(j,k)=dfabs(j,k)-flux(j,1)
            end do

            ! 3.3  Absorption and reflection in the troposphere
            do k=3,kx
                do j=1,ngp
                    tau2(j,k,3)=flux(j,1)*tau2(j,k,3)
                    flux (j,1)=flux(j,1)-tau2(j,k,3)
                    dfabs(j,k)=flux(j,1)
                    flux (j,1)=tau2(j,k,1)*flux(j,1)
                    dfabs(j,k)=dfabs(j,k)-flux(j,1)
                end do
            end do

            do k=2,kx
                do j=1,ngp
                  dfabs(j,k)=dfabs(j,k)+flux(j,2)
                  flux (j,2)=tau2(j,k,2)*flux(j,2)
                  dfabs(j,k)=dfabs(j,k)-flux(j,2)
                end do
            end do

            ! 4. Shortwave upward flux
            ! 4.1  Absorption and reflection at the surface
            do j=1,ngp
                fsfcd(j)  = flux(j,1)+flux(j,2)
                flux(j,1) = flux(j,1)*albsfc(j)
                fsfc(j)   = fsfcd(j)-flux(j,1)
            end do

            ! 4.2  Absorption of upward flux
            do k=kx,1,-1
                do j=1,ngp
                    dfabs(j,k)=dfabs(j,k)+flux(j,1)
                    flux (j,1)=tau2(j,k,1)*flux(j,1)
                    dfabs(j,k)=dfabs(j,k)-flux(j,1)
                    flux (j,1)=flux(j,1)+tau2(j,k,3)
                end do
            end do

            ! 4.3  Net solar radiation = incoming - outgoing
            ftop = ftop - flux(:,1)

            ! 5.  Initialization of longwave radiation model
            ! 5.1  Longwave transmissivity:
            ! function of layer mass, abs. humidity and cloud cover.

            ! Cloud-free levels (stratosphere + PBL)
            k=1
            do j=1,ngp
                deltap=psa(j)*dsig(k)
                tau2(j,k,1)=exp(-deltap*ablwin)
                tau2(j,k,2)=exp(-deltap*ablco2)
                tau2(j,k,3)=1.0_dp
                tau2(j,k,4)=1.0_dp
            end do

            do k=2,kx,kx-2
                do j=1,ngp
                    deltap=psa(j)*dsig(k)
                    tau2(j,k,1)=exp(-deltap*ablwin)
                    tau2(j,k,2)=exp(-deltap*ablco2)
                    tau2(j,k,3)=exp(-deltap*ablwv1*qa(j,k))
                    tau2(j,k,4)=exp(-deltap*ablwv2*qa(j,k))
                end do
            end do

            ! Cloudy layers (free troposphere)
            acloud = cloudc * ablcl2

            do k=3,nl1
               do j=1,ngp
                 deltap=psa(j)*dsig(k)
                 if (k.lt.icltop(j)) then
                   acloud1=acloud(j)
                 else
                   acloud1=ablcl1*cloudc(j)
                 endif
                 tau2(j,k,1)=exp(-deltap*(ablwin+acloud1))
                 tau2(j,k,2)=exp(-deltap*ablco2)
                 tau2(j,k,3)=exp(-deltap*max(ablwv1*qa(j,k),acloud(j)))
                 tau2(j,k,4)=exp(-deltap*max(ablwv2*qa(j,k),acloud(j)))
               end do
            end do

            ! 5.2  Stratospheric correction terms
            eps1=epslw/(dsig(1)+dsig(2))
            do j=1,ngp
                stratc(j,1)=stratz(j)*psa(j)
                stratc(j,2)=eps1*psa(j)
            end do
        end subroutine radsw
end module phy_radsw
