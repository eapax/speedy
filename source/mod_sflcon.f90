module mod_sflcon
    use mod_atparam
    use mod_prec

    implicit none

    private
    public fwind0, ftemp0, fhum0, cdl, cds, chl, chs, vgust, ctday, dtheta,&
        & fstab, hdrag, fhdrag, clambda, clambsn, forog

    !  Constants for surface fluxes
    ! Ratio of near-sfc wind to lowest-level wind
    real(dp) :: fwind0 = 0.95

    ! Weight for near-sfc temperature extrapolation (0-1) :
    !          1 : linear extrapolation from two lowest levels
    !          0 : constant potential temperature ( = lowest level)
    real(dp) :: ftemp0 = 1.0

    ! Weight for near-sfc specific humidity extrapolation (0-1) :
    !            1 : extrap. with constant relative hum. ( = lowest level)
    !            0 : constant specific hum. ( = lowest level)
    real(dp) :: fhum0 = 0.0

    ! Drag coefficient for momentum over land
    real(dp) :: cdl = 2.4e-3

    ! Drag coefficient for momentum over sea
    real(dp) :: cds = 1.0e-3

    ! Heat exchange coefficient over land
    real(dp) :: chl = 1.2e-3

    ! Heat exchange coefficient over sea
    real(dp) :: chs = 0.9e-3

    ! Wind speed for sub-grid-scale gusts
    real(dp) :: vgust = 5.0

    ! Daily-cycle correction (dTskin/dSSRad)
    real(dp) :: ctday = 1.0e-2

    ! Potential temp. gradient for stability correction
    real(dp) :: dtheta = 3.0

    ! Amplitude of stability correction (fraction)
    real(dp) :: fstab = 0.67

    ! Height scale for orographic correction
    real(dp) :: hdrag = 2000.0

    ! Amplitude of orographic correction (fraction)
    real(dp) :: fhdrag = 0.5

    ! Heat conductivity in skin-to-root soil layer
    real(dp) :: clambda = 7.0

    ! Heat conductivity in soil for snow cover = 1
    real(dp) :: clambsn = 7.0

    ! Time-invariant fields (initial. in SFLSET)
    real(dp) :: forog(ix*il)
end module
