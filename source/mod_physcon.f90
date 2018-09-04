module mod_physcon
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Physical constants
    ! Reference pressure
    real(dp), parameter :: p0 = 1.d+5

    ! Gravity accel.
    real(dp), parameter :: gg = 9.81_dp

    ! Gas constant for dry air
    real(dp), parameter :: rd = 287.0_dp

    ! Specific heat at constant pressure
    real(dp), parameter :: cp = 1004.0_dp

    ! Latent heat of condensation, in J/g for consistency with spec.hum. in g/Kg
    real(dp), parameter :: alhc = 2501.0_dp

    ! Latent heat of sublimation
    real(dp), parameter :: alhs = 2801.0_dp

    ! Stefan-Boltzmann constant
    real(dp), parameter :: sbc = 5.67d-8

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
    real(dp), dimension(:), allocatable :: &
            sig, sigl, dsig, pout, grdsig, grdscp
    real(dp), allocatable :: wvi(:,:), sigh(:)
    real(dp), dimension(:), allocatable :: slat, clat

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
end module
