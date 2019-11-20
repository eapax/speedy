module phy_radlw
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_lw_radiation, ini_radlw, truncate_radlw,&
            radlw_transmissivity, radlw_down, radlw_up
    public epslw, emisfc

    ! Number of radiation bands with tau < 1
    integer, parameter :: nband=4

    ! Variables loaded in by namelist
    namelist /lw_radiation/ epslw, emisfc, &
            ablwin, ablwv1, ablwv2, &
            ablcl1, ablcl2

    ! epslw  = fraction of blackbody spectrum absorbed/emitted by PBL only
    type(rpe_var) :: epslw
    ! emisfc = longwave surface emissivity
    type(rpe_var) :: emisfc

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

    ! Time-invariant fields (initial. in radset)
    ! fband  = energy fraction emitted in each LW band = f(T)
    type(rpe_var) :: fband(100:400,4)

    ! Transmissivity and blackbody rad. (updated in radlw)
    ! st4a   = blackbody emission from full and half atmospheric levels
    type(rpe_var), allocatable :: st4a(:,:,:)
    ! flux   = radiative flux in different spectral bands
    type(rpe_var), allocatable :: flux(:,:)
    !  dfabs  = flux of lw rad. absorbed by each atm. layer (3-dim)
    type(rpe_var), allocatable :: dfabs(:,:)
    ! Transmissivity and blackbody radiation
    ! tau2   = transmissivity of atmospheric layers
    ! stratc = stratospheric correction term
    type(rpe_var), allocatable :: tau2(:,:,:), stratc(:,:)

    ! Local derived variables
    type(rpe_var) :: refsfc, eps1

    ! Local copies of mod_physcon variables
    type(rpe_var) :: sbc_1_3, sbc_1_4

    contains
        subroutine setup_lw_radiation(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, lw_radiation)
            write(*, lw_radiation)

            allocate(st4a(ngp,kx,2))
            allocate(flux(ngp,4))
            allocate(dfabs(ngp,kx))
            allocate(tau2(ngp,kx,4))
            allocate(stratc(ngp,2))
        end subroutine setup_lw_radiation

        subroutine ini_radlw()
            ! Calculate local variables for long-wave radiation scheme
            use mod_physcon, only: sbc, dsig

            ! Derived variables
            refsfc=1.0_dp-emisfc
            eps1=epslw/(dsig(1)+dsig(2))
            call radset()

            ! Local copies of mod_physcon
            sbc_1_3 = sbc**(1.0_dp/3.0_dp)
            sbc_1_4 = sbc**(1.0_dp/4.0_dp)
        end subroutine ini_radlw

        subroutine truncate_radlw()
            ! Truncate local variables for long-wave radiation scheme
            ! Namelist variables
            call apply_truncation(epslw)
            call apply_truncation(emisfc)
            call apply_truncation(ablwin)
            call apply_truncation(ablwv1)
            call apply_truncation(ablwv2)
            call apply_truncation(ablcl1)
            call apply_truncation(ablcl2)

            ! Derived variables
            call apply_truncation(eps1)
            call apply_truncation(fband)
            call apply_truncation(refsfc)

            ! Local copies of mod_physcon
            call apply_truncation(sbc_1_3)
            call apply_truncation(sbc_1_4)
        end subroutine truncate_radlw

        subroutine radset()
            ! subroutine radset
            !
            ! Purpose: compute energy fractions in LW bands
            !          as a function of temperature
            integer :: jb, jtemp
            real(dp) :: eps1

            eps1=1.0_dp-epslw

            do jtemp=200,320
                fband(jtemp,2)=(0.148_dp-3.0d-6*(jtemp-247)**2)*eps1
                fband(jtemp,3)=(0.356_dp-5.2d-6*(jtemp-282)**2)*eps1
                fband(jtemp,4)=(0.314_dp+1.0d-5*(jtemp-315)**2)*eps1
                fband(jtemp,1)=eps1 - &
                        (fband(jtemp,2) + fband(jtemp,3) + fband(jtemp,4))
            end do

            do jb=1,4
                do jtemp=100,199
                    fband(jtemp,jb)=fband(200,jb)
                end do
                do jtemp=321,400
                    fband(jtemp,jb)=fband(320,jb)
                end do
            end do
        end subroutine radset

        subroutine radlw_transmissivity(psa, qa, icltop, cloudc)
            ! 5.  Initialization of longwave radiation model

            ! The following variables are derived once per day in other
            ! subroutines and are only used here.
            use mod_fordate, only: ablco2
            use mod_solar, only: stratz
            use mod_physcon, only: dsig

            !  input:   psa    = norm. surface pressure [p/p0] (2-dim)
            type(rpe_var), intent(in) :: psa(ngp)
            !           qa     = specific humidity [g/kg]                (3-dim)
            type(rpe_var), intent(in) :: qa(ngp,kx)
            !           icltop = cloud top level                         (2-dim)
            integer, intent(in) :: icltop(ngp)
            !           cloudc = total cloud cover                       (2-dim)
            type(rpe_var), intent(in) :: cloudc(ngp)

            ! Local variables
            integer :: j, k
            type(rpe_var) :: acloud(ngp), acloud1, deltap

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

            do k=3,kxm
               do j=1,ngp
                 deltap=psa(j)*dsig(k)
                 if (k<icltop(j)) then
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
            do j=1,ngp
                stratc(j,1)=stratz(j)*psa(j)
                stratc(j,2)=eps1*psa(j)
            end do

        end subroutine radlw_transmissivity

        subroutine radlw_down(ta_in,fsfcd)
            !  Purpose: Compute the absorption of longwave radiation
            !           downward flux only
            use mod_physcon, only: wvi

            !  input:  ta     = absolute temperature (3-dim)
            type(rpe_var), intent(in) :: ta_in(ngp,kx)

            !  output:  fsfcd  = downward flux of lw rad. at the sfc.
            type(rpe_var), intent(out) :: fsfcd(ngp)

            ! Local copies of input variables
            type(rpe_var) :: ta(ngp,kx)

            ! Local variables
            integer :: j, jb, k
            type(rpe_var) :: anis, anish, brad, corlw, emis, eps1
            type(rpe_var) :: st3a

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            ta = ta_in + 273.0_dp

            ! 1. Blackbody emission from atmospheric levels.
            ! The linearized gradient of the blakbody emission is computed
            ! from temperatures at layer boundaries, which are interpolated
            ! assuming a linear dependence of T on log_sigma.
            ! Above the first (top) level, the atmosphere is assumed isothermal.

            ! Temperature at level boundaries
            do k=1,kxm
                do j=1,ngp
                    st4a(j,k,1)=ta(j,k)+wvi(k,2)*(ta(j,k+1)-ta(j,k))
                end do
            end do

            ! Mean temperature in stratospheric layers
            do j=1,ngp
                st4a(j,1,2)=rpe_literal(0.75_dp)*ta(j,1) + &
                        rpe_literal(0.25_dp)*st4a(j,1,1)
                st4a(j,2,2)=rpe_literal(0.50_dp)*ta(j,2) + &
                        rpe_literal(0.25_dp)*(st4a(j,1,1)+st4a(j,2,1))
            end do

            ! Temperature gradient in tropospheric layers
            anis =1.0_dp
            anish=rpe_literal(0.5_dp)*anis

            do k=3,kxm
                do j=1,ngp
                    st4a(j,k,2)=anish* &
                            max(st4a(j,k,1)-st4a(j,k-1,1), rpe_literal(0.0_dp))
                end do
            end do

            do j=1,ngp
                st4a(j,kx,2)=anis* &
                        max(ta(j,kx) - st4a(j,kxm,1), rpe_literal(0.0_dp))
            end do

            ! Blackbody emission in the stratosphere
            do k=1,2
                do j=1,ngp
                    st4a(j,k,1)=(sbc_1_4 * st4a(j,k,2))**4
                    st4a(j,k,2)=0.0_dp
                end do
            end do

            ! Blackbody emission in the troposphere
            do k=3,kx
                do j=1,ngp
                    st3a=(sbc_1_3 * ta(j,k))**3
                    st4a(j,k,1)=st3a*ta(j,k)
                    st4a(j,k,2)=rpe_literal(4.0_dp)*st3a*st4a(j,k,2)
                end do
            end do

            ! 2. Initialization of fluxes
            fsfcd = 0.0_dp
            dfabs = 0.0_dp

            ! 3. Emission ad absorption of longwave downward flux.
            !    For downward emission, a correction term depending on the
            !    local temperature gradient and on the layer transmissivity is
            !    added to the average (full-level) emission of each layer.

            ! 3.1  Stratosphere
            k=1
            do jb=1,2
                do j=1,ngp
                    emis=rpe_literal(1.0_dp)-tau2(j,k,jb)
                    brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)+emis*st4a(j,k,2))
                    flux(j,jb)=emis*brad
                    dfabs(j,k)=dfabs(j,k)-flux(j,jb)
                end do
            end do

            flux(:,3:nband) = 0.0_dp

            ! 3.2  Troposphere
            do jb=1,nband
                do k=2,kx
                    do j=1,ngp
                        emis=rpe_literal(1.0_dp)-tau2(j,k,jb)
                        brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1) + &
                                emis*st4a(j,k,2))
                        dfabs(j,k)=dfabs(j,k)+flux(j,jb)
                        flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
                        dfabs(j,k)=dfabs(j,k)-flux(j,jb)
                    end do
                end do
            end do

            ! 3.3 Surface downward flux
            do jb=1,nband
                do j=1,ngp
                    fsfcd(j)=fsfcd(j)+emisfc*flux(j,jb)
                end do
            end do

            ! 3.4 Correction for "black" band (incl. surface reflection)
            eps1=epslw*emisfc
            do j=1,ngp
                corlw=eps1*st4a(j,kx,1)
                dfabs(j,kx)=dfabs(j,kx)-corlw
                fsfcd(j)     =fsfcd(j)     +corlw
            end do

        end subroutine radlw_down

        subroutine radlw_up(ta_in, ts, fsfcd, fsfcu, flx2tend, &
                fsfc, ftop, tt_rlw)
            !  Purpose: Compute the absorption of longwave radiation
            !           upward flux only
            use mod_physcon, only: dsig

            !           ta     = absolute temperature (3-dim)
            type(rpe_var), intent(in) :: ta_in(ngp,kx)
            !           ts     = surface temperature
            type(rpe_var), intent(in) :: ts(ngp)
            !           fsfcd  = downward flux of lw rad. at the sfc.
            type(rpe_var), intent(in) :: fsfcd(ngp)
            !           fsfcu  = surface blackbody emission (upward)
            type(rpe_var), intent(in) :: fsfcu(ngp)
            !           flx2tend = Conversion from fluxes to temperature tendencies
            type(rpe_var), intent(in) :: flx2tend(ngp,kx)

            !  Output:  fsfc   = net upw. flux of lw rad. at the sfc.
            type(rpe_var), intent(out) :: fsfc(ngp)
            !           ftop   = outgoing flux of lw rad. at the top
            type(rpe_var), intent(out) :: ftop(ngp)
            !           tt_rlw = Temperature tendency due to LW radiation
            type(rpe_var), intent(out) :: tt_rlw(ngp,kx)

            ! Local copies of input variables
            type(rpe_var) :: ta(ngp,kx)

            ! Local variables
            integer :: j, jb, k
            type(rpe_var) :: brad, corlw1, corlw2, emis

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            ta = ta_in + 273.0_dp

            fsfc = fsfcu - fsfcd

            do jb=1,nband
                do j=1,ngp
                    flux(j,jb)=fband(nint(ts(j)),jb)*fsfcu(j)+refsfc*flux(j,jb)
                end do
            end do

            ! 4.2  Troposphere

            ! Correction for "black" band
            do j=1,ngp
                dfabs(j,kx)=dfabs(j,kx)+epslw*fsfcu(j)
            end do

            do jb=1,nband
                do k=kx,2,-1
                    do j=1,ngp
                        emis=rpe_literal(1.0_dp)-tau2(j,k,jb)
                        brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1) - &
                                emis*st4a(j,k,2))
                        dfabs(j,k)=dfabs(j,k)+flux(j,jb)
                        flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
                        dfabs(j,k)=dfabs(j,k)-flux(j,jb)
                    end do
                end do
            end do

            ! 4.3  Stratosphere
            k=1
            do jb=1,2
                do j=1,ngp
                    emis=rpe_literal(1.0_dp)-tau2(j,k,jb)
                    brad=fband(nint(ta(j,k)),jb)*(st4a(j,k,1)-emis*st4a(j,k,2))
                    dfabs(j,k)=dfabs(j,k)+flux(j,jb)
                    flux(j,jb)=tau2(j,k,jb)*flux(j,jb)+emis*brad
                    dfabs(j,k)=dfabs(j,k)-flux(j,jb)
                end do
            end do

            ! Correction for "black" band and polar night cooling
            do j=1,ngp
                corlw1=dsig(1)*stratc(j,2)*st4a(j,1,1)+stratc(j,1)
                corlw2=dsig(2)*stratc(j,2)*st4a(j,2,1)
                dfabs(j,1)=dfabs(j,1)-corlw1
                dfabs(j,2)=dfabs(j,2)-corlw2
                ftop(j)   =corlw1+corlw2
            end do

            ! 4.4  Outgoing longwave radiation
            do jb=1,nband
                do j=1,ngp
                    ftop(j)=ftop(j)+flux(j,jb)
                end do
            end do

            ! Convert fluxes to temperature tendencies
            tt_rlw = dfabs*flx2tend
        end subroutine radlw_up
end module phy_radlw
