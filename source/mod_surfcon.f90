module mod_surfcon
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    namelist /surface/ swcap, swwil, sd2sc

    ! Soil moisture parameters
    ! Soil wetness at field capacity (volume fraction)
    real(dp) :: swcap

    ! Soil wetness at wilting point  (volume fraction)
    real(dp) :: swwil

    ! Snow depth (mm water) corresponding to snow cover = 1
    real(dp) :: sd2sc

    ! Land-sea masks (initial. in INBCON)
    ! Original (fractional) land-sea mask
    real(dp), allocatable :: fmask(:,:)
    ! Model-defined land fraction
    real(dp), allocatable :: fmask1(:,:)

    ! Time invariant surface fields
    ! (initial. in INBCON, phis0 initial. in INVARS)
    ! Unfiltered surface geopotential
    real(dp), allocatable :: phi0(:,:)

    ! Spectrally-filtered sfc. geopotential
    real(dp), allocatable :: phis0(:,:)

    ! Bare-land annual-mean albedo
    real(dp), allocatable :: alb0(:,:)

    contains
        subroutine setup_surface(fid)
            integer, intent(in) :: fid
            allocate(fmask(ix,il))
            allocate(fmask1(ix,il))
            allocate(phi0(ix,il))
            allocate(phis0(ix,il))
            allocate(alb0(ix,il))

            read(fid, surface)

            write(*, surface)
        end subroutine setup_surface
end module
