module mod_physcon
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_physcon, init_physcon, truncate_physcon, &
            p0, gg, rd, cp, alhc, alhs, &
            sig, sigl, sigh, dsig, pout, grdsig, grdscp, wvi, slat, clat, &
            sbc_1_3, sbc_1_4

    ! Physical constants
    ! Reference pressure
    real(dp), parameter :: p0 = 1.d+5

    ! Gravity accel.
    real(dp), parameter :: gg_ = 9.81_dp

    ! Gas constant for dry air
    real(dp), parameter :: rd_ = 287.0_dp

    ! Specific heat at constant pressure
    real(dp), parameter :: cp_ = 1004.0_dp

    ! Latent heat of condensation, in J/g for consistency with spec.hum. in g/Kg
    real(dp), parameter :: alhc_ = 2501.0_dp

    ! Latent heat of sublimation
    real(dp), parameter :: alhs_ = 2801.0_dp

    ! Stefan-Boltzmann constant
    real(dp), parameter :: sbc = 5.67d-8

    ! Reduced precision versions
    type(rpe_var) :: gg
    type(rpe_var) :: rd
    type(rpe_var) :: cp
    type(rpe_var) :: alhc
    type(rpe_var) :: alhs

    !   Functions of sigma and latitude (initial. in INPHYS)
    !    sig    = full-level sigma
    !    sigl   = logarithm of full-level sigma
    !    sigh   = half-level sigma
    !    dsig   = layer depth in sigma
    !    pout   = norm. pressure level [p/p0] for post-processing
    !    grdsig = g/(d_sigma p0) : to convert fluxes of u,v,q into d(u,v,q)/dt
    !    grdscp = g/(d_sigma p0 c_p): to convert energy fluxes into dT/dt
    !        Note that grdsig and grdscp are also both multiplied by 3600
    !        (1 hour in seconds) resulting in tendencies being calculated as per
    !        hour when converting from fluxes.
    !    wvi    = weights for vertical interpolation
    !    slat   = sin(lat)
    !    clat   = cos(lat)
    type(rpe_var), dimension(:), allocatable :: &
            sig, sigl, dsig, pout, grdsig, grdscp
    type(rpe_var), allocatable :: wvi(:,:), sigh(:)
    type(rpe_var), dimension(:), allocatable :: slat, clat

    ! Powers of Stefan-Boltzmann constant
    type(rpe_var) :: sbc_1_3, sbc_1_4

    contains
        subroutine setup_physcon()
            allocate(sig(kx))
            allocate(sigl(kx))
            allocate(dsig(kx))
            allocate(pout(kx))
            allocate(grdsig(kx))
            allocate(grdscp(kx))
            allocate(wvi(kx,2))
            allocate(sigh(0:kx))
            allocate(slat(il))
            allocate(clat(il))
        end subroutine setup_physcon

        subroutine init_physcon()
            gg = gg_
            rd = rd_
            cp = cp_
            alhc = alhc_
            alhs = alhs_
            sbc_1_3 = sbc ** (1.0_dp/3.0_dp)
            sbc_1_4 = sbc ** (1.0_dp/4.0_dp)
        end subroutine init_physcon

        subroutine truncate_physcon()
            call apply_truncation(gg)
            call apply_truncation(rd)
            call apply_truncation(cp)
            call apply_truncation(alhc)
            call apply_truncation(alhs)

            call apply_truncation(sig)
            call apply_truncation(sigl)
            call apply_truncation(dsig)
            call apply_truncation(pout)
            call apply_truncation(grdsig)
            call apply_truncation(grdscp)

            call apply_truncation(wvi)
            call apply_truncation(sigh)
            call apply_truncation(slat)
            call apply_truncation(clat)

            call apply_truncation(sbc_1_3)
            call apply_truncation(sbc_1_4)
        end subroutine truncate_physcon
end module
