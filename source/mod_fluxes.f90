module mod_fluxes
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Net heat flux into land sfc.end module
    real(dp), allocatable :: hflux_l(:)

    ! Net heat flux into sea sfc.
    real(dp), allocatable :: hflux_s(:)

    ! Net heat flux into sea-ice sfc.
    real(dp), allocatable :: hflux_i(:)

    contains
        subroutine setup_fluxes()
            allocate(hflux_l(ngp))
            allocate(hflux_s(ngp))
            allocate(hflux_i(ngp))
        end subroutine setup_fluxes

        subroutine ini_fluxes()
            ! Purpose: Initialise fluxes for sea/land models
            !          Set to zero at the start of a day
            !          If initialised part way through a day load from flux file
            use mod_date, only: ihour

            !integer :: j

            if (ihour/=0) then
                !TODO - Implement reading and writing to a fluxes file to start
                !       midway through a day

                ! Read from flux file
                !open (100,file='fluxes.grd',form='unformatted',access='direct',&
                !        recl=8*ngp)
                !read (100,rec=4) (hflux_l(j),j=1,ngp)
                !read (100,rec=14) (hflux_s(j),j=1,ngp)
                !read (100,rec=15) (hflux_i(j),j=1,ngp)
                !close (100)
                stop 'Starting at arbitrary time of day not implemented yet'
            else
                ! Set all daily-mean arrays to zero
                hflux_l(:) = 0.0_dp
                hflux_s(:) = 0.0_dp
                hflux_i(:) = 0.0_dp
            end if
        end subroutine ini_fluxes

        subroutine increment_fluxes()
            ! Purpose: Add up fluxes to provide daily averages
            !          used in sea/land models

            use mod_tsteps, only: nsteps
            use mod_physcon, only: alhc, sbc
            use mod_surfcon, only: fmask1
            use mod_physvar, only: hfluxn, ssrd, shf, evap
            use mod_var_sea, only: tice_am, sice_am
            use mod_fordate, only: albsea, albice
            use phy_radlw, only: emisfc

            real(dp) :: difice(ngp)
            real(dp) :: fland(ngp), esbc, rsteps, sstfr, sstfr4

            fland = reshape(fmask1,(/ngp/))
            rsteps = 1.0_dp/nsteps

            ! SST at freezing point
            sstfr  = 273.2_dp-1.8_dp

            sstfr4 = sstfr**4
            esbc   = emisfc%val*sbc

            ! 2. Store fluxes over land (SI units, all heat fluxes downw.)
            hflux_l(:) = hflux_l(:) + hfluxn(:,1)*rsteps

            ! 3. Store fluxes over sea (SI units, all heat fluxes downw.)
            ! Difference in net (downw.) heat flux between ice and sea surface
            difice(:) = (albsea-albice)*ssrd(:)+ esbc*(sstfr4-tice_am(:)**4)&
                & + shf(:,2)+evap(:,2)*alhc

            hflux_s(:) = hflux_s(:) + rsteps* hfluxn(:,2)
            hflux_i(:) = hflux_i(:) + &
                    rsteps*(hfluxn(:,2)+difice(:)*(1.0_dp-sice_am(:)))

            ! 4.1 Store fluxes for daily-mean output

            ! Multiply net heat fluxes by land or sea fractions
            hfluxn(:,1) = hfluxn(:,1)*fland(:)
            hfluxn(:,2) = hfluxn(:,2)*(1.0_dp-fland(:))

            ! End of flux increment
        end subroutine increment_fluxes

end module mod_fluxes
