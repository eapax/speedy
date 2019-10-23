module phy_radsw
    ! Put this in the documentation please
    use mod_atparam
    use rp_emulator
    use mod_prec

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
            albcl, albcls

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

    ! albcl  = cloud albedo (for cloud cover = 1)
    type(rpe_var) :: albcl
    ! albcls = stratiform cloud albedo (for st. cloud cover = 1)
    type(rpe_var) :: albcls

    ! Local derived variables
    type(rpe_var) :: fband1, fband2
    type(rpe_var), allocatable :: abs1(:), dsig_sw(:)

    contains
        subroutine setup_sw_radiation(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, sw_radiation)
            write(*, sw_radiation)

            allocate(abs1(2:kx))
            allocate(dsig_sw(kx))
        end subroutine setup_sw_radiation

        subroutine ini_radsw()
            ! Calculate local variables for short-wave radiation scheme
            use mod_physcon, only: sig, dsig
            use phy_radlw, only: epslw

            fband2 = 0.05_dp
            fband1 = 1.0_dp-fband2
            abs1(2:kx)=absdry+absaer*sig(2:kx)**2
            dsig_sw(1:kx) = dsig(1:kx)
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
            call apply_truncation(albcl)
            call apply_truncation(albcls)

            ! Locally derived variables
            call apply_truncation(fband1)
            call apply_truncation(fband2)
            call apply_truncation(abs1)

            ! Local copies of mod_physcon
            call apply_truncation(dsig_sw)
        end subroutine truncate_radsw

        subroutine radsw(&
                psa_in, qa_in, icltop, cloudc_in, clstr_in, flx2tend_in, &
                fsfcd, fsfc, ftop, dfabs)
            ! Compute the absorption of shortwave radiation

            ! The following variables are derived once per day in other
            ! subroutines and are only used here. Therefore they are truncated
            ! after being calculated to the precision of radsw
            use mod_fordate, only: albsfc
            use mod_solar, only: fsol, ozone, ozupp, zenit

            !  input:   psa    = norm. surface pressure [p/p0] (2-dim)
            type(rpe_var), intent(in) :: psa_in(ngp)
            !           qa     = specific humidity [g/kg]                (3-dim)
            type(rpe_var), intent(in) :: qa_in(ngp,kx)
            !           icltop = cloud top level                         (2-dim)
            integer, intent(in) :: icltop(ngp)
            !           cloudc = total cloud cover                       (2-dim)
            type(rpe_var), intent(in) :: cloudc_in(ngp)
            !           clstr  = stratiform cloud cover                  (2-dim)
            type(rpe_var), intent(in) :: clstr_in(ngp)
            !         flx2tend = Conversion factor between fluxes and T tendency
            type(rpe_var), intent(in) :: flx2tend_in(ngp,kx)
            !  output:  fsfcd  = downward-only flux of sw rad. at the surface (2-dim)
            type(rpe_var), intent(out) :: fsfcd(ngp)
            !           fsfc   = net (downw.) flux of sw rad. at the surface  (2-dim)
            type(rpe_var), intent(out) :: fsfc(ngp)
            !           ftop   = net (downw.) flux of sw rad. at the atm. top (2-dim)
            type(rpe_var), intent(out) :: ftop(ngp)
            !           dfabs  = flux of sw rad. absorbed by each atm. layer  (3-dim)
            type(rpe_var), intent(out) :: dfabs(ngp,kx)

            ! Local copies of input variables
            type(rpe_var) :: psa(ngp), qa(ngp,kx), cloudc(ngp), clstr(ngp), &
                    flx2tend(ngp,kx)
            ! flux   = radiative flux in different spectral bands
            type(rpe_var) :: flux(ngp,2)

            type(rpe_var) :: tau2(ngp,kx,3)

            ! Local variables
            integer :: j, k
            type(rpe_var) :: acloud(ngp), psaz(ngp), deltap

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            psa = psa_in
            qa = qa_in
            cloudc = cloudc_in
            clstr = clstr_in
            flx2tend = flx2tend_in

            ! 1.  Initialization
            tau2 = 0.0_dp

            do j=1,ngp
                !fk-- change to ensure only icltop <= kx used
                if(icltop(j)<=kx) then
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
                acloud(j)=cloudc(j)*min(abscl1*qa(j,kxm),abscl2)
            end do

            do j=1,ngp
                deltap=psaz(j)*dsig_sw(1)
                tau2(j,1,1)=exp(-deltap*absdry)
            end do

            do k=2,kxm

                do j=1,ngp
                    deltap=psaz(j)*dsig_sw(k)
                    if (k>=icltop(j)) then
                        tau2(j,k,1)=exp(-deltap*(abs1(k)+abswv1*qa(j,k)+acloud(j)))
                    else
                      tau2(j,k,1)=exp(-deltap*(abs1(k)+abswv1*qa(j,k)))
                    endif
                end do
            end do


            do j=1,ngp
                deltap=psaz(j)*dsig_sw(kx)
                tau2(j,kx,1)=exp(-deltap*(abs1(kx)+abswv1*qa(j,kx)))
            end do

            do k=2,kx
                do j=1,ngp
                  deltap=psaz(j)*dsig_sw(k)
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

            ! Convert SW fluxes to temperature tendencies
            dfabs = dfabs*flx2tend
        end subroutine radsw
end module phy_radsw
