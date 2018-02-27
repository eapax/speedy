module mod_surfcon
    use mod_atparam

    implicit none

    namelist /surface/ swcap, swwil, sd2sc

    ! Land-sea masks (initial. in INBCON)
    ! Original (fractional) land-sea mask
    real, allocatable :: fmask(:,:)
    ! Model-defined land fraction
    real, allocatable :: fmask1(:,:)
									
    ! Time invariant surface fields 
    ! (initial. in INBCON, phis0 initial. in INVARS)
    ! Unfiltered surface geopotential
    real, allocatable :: phi0(:,:)

    ! Spectrally-filtered sfc. geopotential
    real, allocatable :: phis0(:,:)

    ! Bare-land annual-mean albedo
    real, allocatable :: alb0(:,:)

    ! Soil moisture parameters
    ! Soil wetness at field capacity (volume fraction)
    real :: swcap = 0.30 

    ! Soil wetness at wilting point  (volume fraction)
    real :: swwil = 0.17

    ! Snow depth (mm water) corresponding to snow cover = 1
    real :: sd2sc = 60.0

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
