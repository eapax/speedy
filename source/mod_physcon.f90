module mod_physcon
    use mod_atparam

    implicit none

    ! Physical constants
    ! Reference pressure
    real, parameter :: p0 = 1.e+5

    ! Gravity accel.
    real, parameter :: gg = 9.81

    ! Gas constant for dry air
    real, parameter :: rd = 287.

    ! Specific heat at constant pressure
    real, parameter :: cp = 1004.

    ! Latent heat of condensation, in J/g for consistency with spec.hum. in g/Kg
    real, parameter :: alhc = 2501.0

    ! Latent heat of sublimation
    real, parameter :: alhs = 2801.0

    ! Stefan-Boltzmann constant
    real, parameter :: sbc = 5.67e-8

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
    real, dimension(:), allocatable :: sig, sigl, dsig, pout, grdsig, grdscp
    real, allocatable :: wvi(:,:), sigh(:)
    real, dimension(:), allocatable :: slat, clat

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
