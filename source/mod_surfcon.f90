module mod_surfcon
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    namelist /surface/ swcap, swwil, sd2sc

    ! Soil moisture parameters
    ! Soil wetness at field capacity (volume fraction)
    type(rpe_var) :: swcap

    ! Soil wetness at wilting point  (volume fraction)
    type(rpe_var) :: swwil

    ! Snow depth (mm water) corresponding to snow cover = 1
    type(rpe_var) :: sd2sc

    ! Land-sea masks (initial. in INBCON)
    ! Original (fractional) land-sea mask
    type(rpe_var), allocatable :: fmask(:,:)
    ! Model-defined land fraction
    type(rpe_var), allocatable :: fmask1(:,:)
									
    ! Time invariant surface fields 
    ! (initial. in INBCON, phis0 initial. in INVARS)
    ! Unfiltered surface geopotential
    type(rpe_var), allocatable :: phi0(:,:)

    ! Spectrally-filtered sfc. geopotential
    type(rpe_var), allocatable :: phis0(:,:)

    ! Bare-land annual-mean albedo
    type(rpe_var), allocatable :: alb0(:,:)

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

        subroutine init_surfcon
            swcap = swcap
            swwil = swwil
            sd2sc = sd2sc
        end subroutine

        subroutine truncate_surfcon()
            fmask = fmask
            fmask1 = fmask1
            phi0 = phi0
            phis0 = phis0
            alb0 = alb0
        end subroutine truncate_surfcon
end module
