module mod_surfcon
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public fmask, fmask1, phi0, phis0, alb0, swcap, swwil, sd2sc

    ! Land-sea masks (initial. in INBCON)
    ! Original (fractional) land-sea mask
    type(rpe_var) :: fmask(ix,il)
    ! Model-defined land fraction
    type(rpe_var) :: fmask1(ix,il)
									
    ! Time invariant surface fields 
    ! (initial. in INBCON, phis0 initial. in INVARS)
    ! Unfiltered surface geopotential
    type(rpe_var) :: phi0(ix,il)

    ! Spectrally-filtered sfc. geopotential
    type(rpe_var) :: phis0(ix,il)

    ! Bare-land annual-mean albedo
    type(rpe_var) :: alb0(ix,il)

    ! Soil moisture parameters
    ! Soil wetness at field capacity (volume fraction)
    real :: swcap = 0.30 

    ! Soil wetness at wilting point  (volume fraction)
    real :: swwil = 0.17

    ! Snow depth (mm water) corresponding to snow cover = 1
    real :: sd2sc = 60.0
end module
