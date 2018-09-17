module phy_radlw
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public setup_lw_radiation, ini_radlw, radlw_down, radlw_up
    public epslw, emisfc

    ! Number of radiation bands with tau < 1
    integer, parameter :: nband=4

    ! Variables loaded in by namelist
    namelist /lw_radiation/ epslw, emisfc

    ! epslw  = fraction of blackbody spectrum absorbed/emitted by PBL only
    real(dp) :: epslw
    ! emisfc = longwave surface emissivity
    real(dp) :: emisfc

    ! Time-invariant fields (initial. in radset)
    ! fband  = energy fraction emitted in each LW band = f(T)
    real(dp) :: fband(100:400,4)

    ! Transmissivity and blackbody rad. (updated in radlw)
    ! st4a   = blackbody emission from full and half atmospheric levels
    real(dp), allocatable :: st4a(:,:,:)
    ! flux   = radiative flux in different spectral bands
    real(dp), allocatable :: flux(:,:)
    !  dfabs  = flux of lw rad. absorbed by each atm. layer (3-dim)
    real(dp), allocatable :: dfabs(:,:)

    ! Local derived variables
    real(dp) :: refsfc

    contains
        subroutine setup_lw_radiation(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, lw_radiation)
            write(*, lw_radiation)

            allocate(st4a(ngp,kx,2))
            allocate(flux(ngp,4))
            allocate(dfabs(ngp,kx))
        end subroutine setup_lw_radiation

        subroutine ini_radlw()
            ! Calculate local variables for long-wave radiation scheme

            ! Derived variables
            refsfc=1.0_dp-emisfc
            call radset()
        end subroutine ini_radlw

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

        subroutine radlw_down(ta, fsfcd)
            !  Purpose: Compute the absorption of longwave radiation
            !           downward flux only
            use mod_physcon, only: sbc, wvi
            use mod_physvar, only: tau2

            !  input:  ta     = absolute temperature (3-dim)
            real(dp), intent(in) :: ta(ngp,kx)

            !  output:  fsfcd  = downward flux of lw rad. at the sfc.
            real(dp), intent(out) :: fsfcd(ngp)

            ! Local variables
            integer :: j, jb, k
            real(dp) :: anis, anish, brad, corlw, emis, &
                    eps1, esbc
            real(dp) :: st3a

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
                st4a(j,1,2)=0.75_dp*ta(j,1) + &
                        0.25_dp*st4a(j,1,1)
                st4a(j,2,2)=0.50_dp*ta(j,2) + &
                        0.25_dp*(st4a(j,1,1)+st4a(j,2,1))
            end do

            ! Temperature gradient in tropospheric layers
            anis =1.0_dp
            anish=0.5_dp*anis

            do k=3,kxm
                do j=1,ngp
                    st4a(j,k,2)=anish* &
                            max(st4a(j,k,1)-st4a(j,k-1,1), 0.0_dp)
                end do
            end do

            do j=1,ngp
                st4a(j,kx,2)=anis* &
                        max(ta(j,kx) - st4a(j,kxm,1), 0.0_dp)
            end do

            ! Blackbody emission in the stratosphere
            do k=1,2
                do j=1,ngp
                    st4a(j,k,1)=sbc*st4a(j,k,2)**4
                    st4a(j,k,2)=0.0_dp
                end do
            end do

            ! Blackbody emission in the troposphere
            do k=3,kx
                do j=1,ngp
                    st3a=sbc*ta(j,k)**3
                    st4a(j,k,1)=st3a*ta(j,k)
                    st4a(j,k,2)=4.0_dp*st3a*st4a(j,k,2)
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
                    emis=1.0_dp-tau2(j,k,jb)
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
                        emis=1.0_dp-tau2(j,k,jb)
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

        subroutine radlw_up(ta, ts, fsfcd, fsfcu, flx2tend, &
                fsfc, ftop, tt_rlw)
            !  Purpose: Compute the absorption of longwave radiation
            !           upward flux only

            ! Variables calculated only on sub set of timesteps in radsw
            use mod_physcon, only: dsig
            use mod_physvar, only: tau2, stratc

            !           ta     = absolute temperature (3-dim)
            real(dp), intent(in) :: ta(ngp,kx)
            !           ts     = surface temperature
            real(dp), intent(in) :: ts(ngp)
            !           fsfcd  = downward flux of lw rad. at the sfc.
            real(dp), intent(in) :: fsfcd(ngp)
            !           fsfcu  = surface blackbody emission (upward)
            real(dp), intent(in) :: fsfcu(ngp)
            !           flx2tend = Conversion from fluxes to temperature tendencies
            real(dp), intent(in) :: flx2tend(ngp,kx)

            !  Output:  fsfc   = net upw. flux of lw rad. at the sfc.
            real(dp), intent(out) :: fsfc(ngp)
            !           ftop   = outgoing flux of lw rad. at the top
            real(dp), intent(out) :: ftop(ngp)
            !           tt_rlw = Temperature tendency due to LW radiation
            real(dp), intent(out) :: tt_rlw(ngp,kx)

            ! Local variables
            integer :: j, jb, k
            real(dp) :: brad, corlw, corlw1, corlw2, emis

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
                        emis=1.0_dp-tau2(j,k,jb)
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
                    emis=1.0_dp-tau2(j,k,jb)
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