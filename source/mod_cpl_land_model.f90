module mod_cpl_land_model
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public truncate_land_model, land_model_init, land_model
    public vland_input, vland_output

    ! 1./heat_capacity (land)
    type(rpe_var) :: rhcapl(ix,il)           

    ! 1./dissip_time (land)
    type(rpe_var) :: cdland(ix,il)           

    ! Input and output land variables exchanged by coupler
    ! Land model input variables
    type(rpe_var) :: vland_input(ix*il,4)            

    ! Land model output variables
    type(rpe_var) :: vland_output(ix*il,2)           

    contains
        subroutine truncate_land_model()
            rhcapl = rhcapl
            cdland = cdland
            vland_input = vland_input
            vland_output = vland_output
        end subroutine

        subroutine land_model_init(fmask_l,alb0) 
            ! subroutine land_model_init (fmask_l,alb0)
            !
            ! purpose : initialization of land model
            ! initialized common blocks: land_mc
            
            ! Input variables
            ! Land mask (fraction of land)
            type(rpe_var), intent(in) :: fmask_l(ix,il)            
            ! Annual-mean albedo
            type(rpe_var), intent(in) :: alb0(ix,il)            
        
            ! Auxiliary variables
            integer :: i, j
            type(rpe_var) :: dmask(ix,il)           ! domain mask
            type(rpe_var) :: depth_soil, depth_lice, tdland, hcapl, hcapli, flandmin
        
            ! 1. Set heat capacities and dissipation times for 
            !    soil and ice-sheet layers 
        
            ! Model parameters (default values)
        
            ! Soil layer depth (m)
            depth_soil = 1.0
        
            ! Land-ice depth (m)
            depth_lice = 5.0
        
            ! Dissipation time (days) for land-surface temp. anomalies
            tdland  = 40.
        
            ! Minimum fraction of land for the definition of anomalies
            flandmin = 1./3.
        
            ! Reset model parameters
            include "cls_inland.h"
        
            ! Heat capacities per m^2 (depth*heat_cap/m^3)
            hcapl  = depth_soil*2.50e+6
            hcapli = depth_lice*1.93e+6
        
            ! 2. Compute constant fields
            ! Set domain mask (blank out sea points)
            dmask(:,:) = 1.
        
            do j=1,il
                do i=1,ix
                    if (fmask_l(i,j).lt.flandmin) dmask(i,j) = 0
                end do
            end do
        
            ! Set time_step/heat_capacity and dissipation fields
            do j=1,il
                do i=1,ix
                    if (alb0(i,j).lt.0.4) then
                        rhcapl(i,j) = 86400./hcapl
                    else
                        rhcapl(i,j) = 86400./hcapli
                    endif
                end do
            end do
        
            cdland(:,:) = dmask(:,:)*tdland/(1.+dmask(:,:)*tdland)
        end
        
        subroutine land_model 
            ! subroutine land_model
            !
            ! purpose : integrate slab land-surface model for one day
        							
        
            ! Input variables:
            type(rpe_var) :: stl0(ix*il)    ! land temp. at initial time
            type(rpe_var) :: hfland(ix*il)    ! land sfc. heat flux between t0 and t1
            type(rpe_var) :: stlcl1(ix*il)    ! clim. land temp. at final time 
        
            ! Output variables
            type(rpe_var) :: stl1(ix*il)     ! land temp. at final time
        
            ! Auxiliary variables
            type(rpe_var) :: hflux(ix*il)   ! net sfc. heat flux
            type(rpe_var) :: tanom(ix*il)   ! sfc. temperature anomaly

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
            tanom = reshape(cdland, (/ ix*il /))*&
                & (tanom+reshape(rhcapl, (/ ix*il /))*hflux)
 
            ! Full SST at final time
            stl1 = tanom + stlcl1

            vland_output(:,1) = stl1
        end
end
