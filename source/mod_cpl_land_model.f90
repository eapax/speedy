module mod_cpl_land_model
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public vland_input, vland_output
    public setup_land, land_model_init, land_model

    namelist /land/ depth_soil, depth_lice, tdland, flandmin

    ! Derived model constants set up in land_model_init
    ! 1./heat_capacity (land)
    real(dp), allocatable :: rhcapl(:, :)
    ! 1./dissip_time (land)
    real(dp), allocatable :: cdland(:, :)

    ! Input and output land variables exchanged by coupler
    ! Land model input variables
    real(dp), allocatable :: vland_input(:, :)
    ! Land model output variables
    real(dp), allocatable :: vland_output(:, :)

    ! Namelist parameters used to set up model constants
    ! Soil layer depth (m)
    real(dp) :: depth_soil
    ! Land-ice depth (m)
    real(dp) :: depth_lice
    ! Dissipation time (days) for land-surface temp. anomalies
    real(dp) :: tdland
    ! Minimum fraction of land for the definition of anomalies (denominator)
    real(dp) :: flandmin

    contains
        subroutine setup_land(fid)
            integer, intent(in) :: fid

            allocate(rhcapl(ix,il))
            allocate(cdland(ix,il))
            allocate(vland_input(ngp,3))
            allocate(vland_output(ngp,1))

            read(fid, land)
            flandmin = 1.0_dp/flandmin

            write(*, land)
        end subroutine setup_land

        subroutine land_model_init(fmask_l,alb0)
            ! subroutine land_model_init (fmask_l,alb0)
            !
            ! purpose : initialization of land model
            ! initialized common blocks: land_mc

            ! Input variables
            ! Land mask (fraction of land)
            real(dp), intent(in) :: fmask_l(ix,il)
            ! Annual-mean albedo
            real(dp), intent(in) :: alb0(ix,il)

            ! Auxiliary variables
            integer :: i, j
            real(dp) :: dmask(ix,il)           ! domain mask
            real(dp) :: hcapl, hcapli

            ! 1. Set heat capacities and dissipation times for
            !    soil and ice-sheet layers

            ! Heat capacities per m^2 (depth*heat_cap/m^3)
            hcapl  = depth_soil*2.50d+6
            hcapli = depth_lice*1.93d+6

            ! 2. Compute constant fields
            ! Set domain mask (blank out sea points)
            dmask(:,:) = 1.0_dp

            do j=1,il
                do i=1,ix
                    if (fmask_l(i,j).lt.flandmin) dmask(i,j) = 0
                end do
            end do

            ! Set time_step/heat_capacity and dissipation fields
            do j=1,il
                do i=1,ix
                    if (alb0(i,j).lt.0.4_dp) then
                        rhcapl(i,j) = 86400.0_dp/hcapl
                    else
                        rhcapl(i,j) = 86400.0_dp/hcapli
                    endif
                end do
            end do

            cdland(:,:) = dmask(:,:)*tdland/(1.0_dp+dmask(:,:)*tdland)
        end subroutine land_model_init

        subroutine land_model()
            ! subroutine land_model
            !
            ! purpose : integrate slab land-surface model for one day


            ! Input variables:
            real(dp) :: stl0(ngp)    ! land temp. at initial time
            real(dp) :: hfland(ngp)    ! land sfc. heat flux between t0 and t1
            real(dp) :: stlcl1(ngp)    ! clim. land temp. at final time

            ! Output variables
            real(dp) :: stl1(ngp)     ! land temp. at final time

            ! Auxiliary variables
            real(dp) :: hflux(ngp)   ! net sfc. heat flux
            real(dp) :: tanom(ngp)   ! sfc. temperature anomaly

            ! Initialise variables
            stl0 = vland_input(:,1)
            hfland = vland_input(:,2)
            stlcl1 = vland_input(:,3)

            ! 1. Land-surface (soil/ice-sheet) layer

            ! Net heat flux
            ! (snow correction to be added?)
            hflux = hfland

            ! Anomaly w.r.t final-time climatological temp.
            tanom = stl0 - stlcl1
            ! Time evoloution of temp. anomaly
            tanom = reshape(cdland, (/ ngp /))*&
                & (tanom+reshape(rhcapl, (/ ngp /))*hflux)

            ! Full SST at final time
            stl1 = tanom + stlcl1

            vland_output(:,1) = stl1
        end subroutine land_model
end module mod_cpl_land_model
