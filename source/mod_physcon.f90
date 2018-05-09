module mod_physcon
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    ! Physical constants
    ! Reference pressure
    real(dp), parameter :: p0_ = 1.e+5

    ! Gravity accel.
    real(dp), parameter :: gg_ = 9.81

    ! Gas constant for dry air
    real(dp), parameter :: rd_ = 287.

    ! Specific heat at constant pressure
    real(dp), parameter :: cp_ = 1004.

    ! Latent heat of condensation, in J/g for consistency with spec.hum. in g/Kg
    real(dp), parameter :: alhc_ = 2501.0

    ! Latent heat of sublimation
    real(dp), parameter :: alhs_ = 2801.0

    ! Stefan-Boltzmann constant
    real(dp), parameter :: sbc_ = 5.67e-8

    ! Reduced precision versions
    type(rpe_var) :: p0
    type(rpe_var) :: gg
    type(rpe_var) :: rd
    type(rpe_var) :: cp
    type(rpe_var) :: alhc
    type(rpe_var) :: alhs
    type(rpe_var) :: sbc

    !   Functions of sigma and latitude (initial. in INPHYS)
    !    sig    = full-level sigma 
    !    sigl   = logarithm of full-level sigma
    !    sigh   = half-level sigma
    !    dsig   = layer depth in sigma
    !    pout   = norm. pressure level [p/p0] for post-processing
    !    grdsig = g/(d_sigma p0) : to convert fluxes of u,v,q into d(u,v,q)/dt
    !    grdscp = g/(d_sigma p0 c_p): to convert energy fluxes into dT/dt
    !    wvi    = weights for vertical interpolation
    !    slat   = sin(lat)
    !    clat   = cos(lat)
    type(rpe_var), dimension(:), allocatable :: &
            sig, sigl, dsig, pout, grdsig, grdscp
    type(rpe_var), allocatable :: wvi(:,:), sigh(:)
    type(rpe_var), dimension(:), allocatable :: slat, clat

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
            p0 = p0_
            gg = gg_
            rd = rd_
            cp = cp_
            alhc = alhc_
            alhs = alhs_
            sbc = sbc_
        end subroutine init_physcon

        subroutine truncate_physcon()
            sig = sig
            sigl = sigl
            dsig = dsig
            pout = pout
            grdsig = grdsig
            grdscp = grdscp
            wvi = wvi
            sigh = sigh
            slat = slat
            clat = clat
        end subroutine truncate_physcon
end module
