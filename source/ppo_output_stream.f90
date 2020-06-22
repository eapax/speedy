module ppo_output_stream

    use mod_atparam
    use mod_dynvar, only: vor, div, t, ps, tr, phi, phis
    use mod_physvar
    use mod_physcon, only: gg, sig, pout
    use mod_cli_sea, only: deglat_s
    use mod_date, only: imonth, month_start
    use humidity, only: zero_C
    use spectral, only: uvspec, grid
    use ppo_plevs, only: pressure_levels, np
    use mod_prec, only: sp, dp

    use netcdf

    implicit none

    private
    public initialise_output, update_output, close_output, check

    ! Generic type used to describe an output stream
    type output_stream
        !
        character (len=100) :: filename

        ! Spectral and Grid variables need to be treated differently in
        ! separate streams due to their different grids and the use of complex
        ! numbers
        logical :: spectral

        ! Flag for interpolating variables to pressure levels
        logical :: plevs

        ! Flags for determining output frequency
        ! nstpinc: How often the output variables are incremented
        ! nstpout: How often the variables are output and increments reset
        ! nstpopen: How often a new file is opened
        integer :: nstpinc, nstpout, nstpopen

        ! The ID used in open/read/write statements
        integer :: file_ID
        integer, allocatable :: nc_var_ID(:)

        ! An array of integers corresponding to the variables being output
        integer :: nvars
        integer, allocatable :: var_ID(:)

        ! The record counter used for direct access in read/write statements
        integer :: rec = 1
        integer :: rec_varid
    end type output_stream

    ! Array containing all output streams. Allocated during model set up
    integer :: nstreams
    type(output_stream), allocatable :: streams(:)

    contains
        ! Initialise output streams from the input .txt file
        subroutine initialise_output()
            namelist /output/ nstreams

            namelist /output_file/ filename, spectral, plevs, nstpinc, &
                    nstpout, nstpopen, nvars
            character (len=100) :: filename
            logical :: spectral, plevs
            integer :: nstpinc, nstpout, nstpopen
            integer :: nvars

            namelist /variables/ var_IDs
            integer, allocatable :: var_IDs(:)

            integer :: n

            ! Read output parameters from input text file
            open(99, file='output_requests.nml')
            read(99, output)
            allocate(streams(nstreams))

            ! Setup output stream for each set of outputs
            do n=1, nstreams
                read(99, output_file)
                allocate(var_IDs(nvars))
                read(99, variables)
                streams(n) = init_output_stream(filename, spectral, plevs, &
                        nstpinc, nstpout, nstpopen, nvars, var_IDs)
                deallocate(var_IDs)
            end do

            ! Output 0'th timestep variables
            call update_output(0)
        end subroutine

        ! Update all of the output streams
        subroutine update_output(istep)
            integer, intent(in) :: istep
            integer :: n, m, p, k
            logical :: l_transform_field(5)

            ! Perform spectral transforms to update gridpoint variables if they
            ! are required
            ! Check which fields need to be updated
            l_transform_field = .false.
            do n=1, nstreams
                ! Check all output streams that need to be updated and require
                ! grid point fields
                if (.not. streams(n)%spectral .and. &
                        xmod(istep, streams(n)%nstpinc)) then
                    ! Check for every prognostic variable in the output stream
                    do m=1, streams(n)%nvars
                        do p=1, 5
                            if (streams(n)%var_ID(m)==p) then
                                l_transform_field(p) = .true.
                            end if
                        end do
                    end do
                end if
            end do

            ! Perform gridpoint transforms
            ! ug1, vg1, tg1, qg1, phig1, pslg1
            ! u and/or v
            if (l_transform_field(1) .or. l_transform_field(2)) then
                do k=1,kx
                    call uvspec(vor(:,:,k,1), div(:,:,k,1), ug1(:,k), vg1(:,k))
                end do

                ug1 = ug1 / 3600.0_dp
                vg1 = vg1 / 3600.0_dp
            end if

            ! Temperature
            if (l_transform_field(3) .or. l_transform_field(5)) then
                do k=1,kx
                    call grid(t(:,:,k,1), tg1(:,k), 1)
                end do

                ! Convert temperature to Kelvin
                tg1 = tg1 + zero_C
            end if

            ! Humidity
            if (l_transform_field(4)) then
                do k=1,kx
                    call grid(tr(:,:,k,1,1), qg1(:,k), 1)
                end do
            end if

            ! Geopotential
            if (l_transform_field(5)) then
                call geop(1)
                do k=1,kx
                    call grid(phi(:,:,k), phig1(:,k), 1)
                end do
            end if

            do n=1, nstreams
                call update_output_stream(streams(n), istep)
            end do
        end subroutine

        subroutine close_output()
            integer :: n

            do n=1, nstreams
                call check( nf90_close(streams(n)%file_ID) )
            end do
        end subroutine close_output

        !
        function init_output_stream(filename, spectral, plevs, nstpinc, &
                nstpout, nstpopen, nvars, var_ID) &
                result(stream)
            character(len=*), intent(in) :: filename
            logical, intent(in) :: spectral, plevs
            integer, intent(in) :: nstpinc, nstpout, nstpopen
            integer, intent(in) :: nvars, var_ID(:)
            type(output_stream) :: stream

            ! Pass parameters to the new object
            stream%filename = filename
            stream%spectral = spectral
            stream%plevs = plevs
            stream%nstpinc = nstpinc
            stream%nstpout = nstpout
            stream%nstpopen = nstpopen
            stream%nvars = nvars
            allocate(stream%nc_var_ID(stream%nvars))
            allocate(stream%var_ID(stream%nvars))
            stream%var_ID = var_ID
        end function

        ! Update called once per timestep for each output stream.
        ! The subroutine determines actions based on output_stream parameters
        ! nstpinc: How often the output variables are increments
        ! nstpout: How often the variables are output and increments reset
        ! nstpopen: How often a new file is opened
        ! Each variable is an integer number of timesteps or, if the integer
        ! is negative, number of months.
        subroutine update_output_stream(stream, istep)
            type(output_stream), intent(inout) :: stream
            integer, intent(in) :: istep

            if (xmod(istep, stream%nstpopen)) call init_nc(stream)

            if (xmod(istep, stream%nstpinc)) call incr_output_stream(stream)

            if (xmod(istep, stream%nstpout)) call write_output_stream(stream, istep)
        end subroutine

        ! Check whether the timestep matches the frequency by using mod in
        ! terms of timesteps if positive and mod in terms of months if negative
        ! If zero returns false
        function xmod(istep, frequency)
            integer, intent(in) :: istep
            integer, intent(in) :: frequency
            logical :: xmod

            if (frequency>0) then
                xmod = (mod(istep, frequency)==0)

            else if (frequency < 0) then
                ! Output on the first day of each month
                xmod = (month_start() .and. mod(imonth, -frequency)==0)
            else
                xmod = .false.
            end if
        end function

        ! todo Increment the output on substeps relative to how often it is
        ! written to file
        subroutine incr_output_stream(stream)
            type(output_stream), intent(inout) :: stream
        end subroutine

        ! Write each of the variables associated with the given stream to file.
        ! Acts a wrapper to separate write functions for spectral and grid
        ! variables.
        subroutine write_output_stream(stream, istep)
            use mod_tsteps, only: delt

            type(output_stream), intent(inout) :: stream
            integer, intent(in) :: istep

            ! Write the next entry in the time dimension
            call check( nf90_put_var(stream%file_ID, stream%rec_varid, (/istep*delt%val/), &
                                     start=(/stream%rec/), count=(/1/)) )

            if (stream%spectral) then
                call write_spectral(stream)
            else
                call write_grid(stream)
            end if

            stream%rec = stream%rec + 1
        end subroutine

        subroutine write_spectral(stream)
            type(output_stream), intent(inout) :: stream
            complex(dp) :: output(mx, nx, kx)
            real(sp) :: re_output(mx, nx, kx), im_output(mx, nx, kx)
            integer :: n, k

            do n=1, stream%nvars
                output = fetch_output_spectral(stream%var_ID(n))

                ! Write real and imaginary parts as separate variables
                re_output = REAL(REAL (output))
                im_output = REAL(AIMAG(output))

                do k=kx, 1, -1
                    write(stream%file_ID, rec=stream%rec) re_output(:, :, k)
                    stream%rec = stream%rec + 1
                end do

                do k=kx, 1, -1
                    write(stream%file_ID, rec=stream%rec) im_output(:, :, k)
                    stream%rec = stream%rec + 1
                end do
            end do
        end subroutine

        subroutine write_grid(stream)
            type(output_stream), intent(inout) :: stream
            real(dp) :: output(ngp, kx)
            real(dp) :: output_p(ngp, np)
            logical :: l_3d
            integer :: n, k

            do n=1, stream%nvars
                call fetch_output_grid(stream%var_ID(n), output, l_3d)

                if (l_3d) then
                    if(stream%plevs) then
                        ! Interpolate output to pressure levels
                        call pressure_levels(stream%var_ID(n), output, output_p)
                        do k=np, 1, -1
                            call check( nf90_put_var(stream%file_ID, stream%nc_var_ID(n), output_p(:, k), &
                                                     start = (/ 1, 1, np-k+1, stream%rec /), &
                                                     count = (/ ix, il, 1, 1/)) )
                        end do
                    else
                        ! Otherwise write model level output
                        if (stream%var_ID(n)==5) then
                            output = output / gg !m
                        end if

                        ! For some reason the height levels need to be written backwards
                        do k=kx, 1, -1
                            call check( nf90_put_var(stream%file_ID, stream%nc_var_ID(n), output(:, k), &
                                                     start = (/ 1, 1, kx-k+1, stream%rec /), &
                                                     count = (/ ix, il, 1, 1/)) )
                        end do
                    end if
                else
                    call check( nf90_put_var(stream%file_ID, stream%nc_var_ID(n), output(:, 1), &
                                             start = (/ 1, 1, stream%rec /), &
                                             count = (/ ix, il, 1/)) )
                end if
            end do
        end subroutine write_grid

        ! Get the variable corresponding to the varID in single bit precision
        ! Essentially a look up for all the variables in mod_dynvar
        ! TODO implement 2d variables in same interface
        function fetch_output_spectral(varID) result(output)
            integer :: varID
            complex(dp) :: output(mx, nx, kx)

            select case(varID)
                ! vorticity = vor(mx, nx, kx, 2)
                case(1)
                output = vor(:, :, :, 2)

                ! divergence = div(mx, nx, kx, 2)
                case(2)
                output = div(:, :, :, 2)

                ! temperature = T(mx, nx, kx, 2)
                case(3)
                output = T(:, :, :, 2)

                ! Log of (norm.) sfc pressure (p_s/p0) = ps(mx, nx, 2)
                !case(4)
                !output(:, :, 1) = ps(:, :, 2)

                ! Tracers = tr(mx, nx, kx, 2, ntr)
                case(5)
                output = tr(:, :, :, 2, 1)

                ! Geopotential = phi(mx, nx, kx)
                case(6)
                output = phi

                ! surface geopotential = phis(mx, nx)
                !case(7)
                !output(:, : , 1) = phis

                case default
                print *, 'Variable no.', varID, ' unavailable for output'
            end select
        end function

        ! Get the variable corresponding to the varID in single bit precision
        ! Essentially a look up for all the variables in mod_physvar
        subroutine fetch_output_grid(varID, output, l_3d)
            use mod_fordate, only: alb_l, alb_s, albsfc, snowc
            use mod_solar, only: fsol, ozone, ozupp, zenit, stratz
            use mod_var_land, only: stl_am, soilw_am
            use mod_var_sea, only: sst_am, ssti_om
            use mod_surfcon, only: phis0, fmask1

            integer :: varID
            real(dp), intent(out) :: output(ngp, kx)

            ! Flag for whether the output is 3d (uses all kx levels)
            logical, intent(out) :: l_3d

            l_3d = .true.

            select case(varID)
                ! ug1    = u-wind
                case(1)
                output = ug1

                ! vg1    = v-wind
                case(2)
                output = vg1

                ! tg1    = abs. temperature
                case(3)
                output = tg1

                ! qg1    = specific humidity (g/kg)
                case(4)
                output = qg1 * 1.0d-3  ! kg/kg

                ! phig1  = geopotential
                case(5)
                output = phig1

                ! pslg1  = log. of surface pressure
                case(6)
                output(:, 1) = pslg1
                l_3d = .false.

                ! se     = dry static energy
                case(7)
                output = se

                ! rh     = relative humidity
                case(8)
                output = rh

                ! qsat   = saturation specific humidity (g/kg)
                case(9)
                output = qsat

                ! psg    = surface pressure
                case(10)
                output(:, 1) = psg
                l_3d = .false.

                ! ts     = surface temperature
                case(11)
                output(:, 1) = ts
                l_3d = .false.

                ! tskin  = skin temperature
                case(12)
                output(:, 1) = tskin
                l_3d = .false.

                ! u0     = near-surface u-wind
                case(13)
                output(:, 1) = u0
                l_3d = .false.

                ! v0     = near-surface v-wind
                case(14)
                output(:, 1) = v0
                l_3d = .false.

                ! t0     = near-surface air temperature
                case(15)
                output(:, 1) = t0
                l_3d = .false.

                ! q0     = near-surface specific humidity (g/kg)
                case(16)
                output(:, 1) = q0
                l_3d = .false.

                ! cloudc = total cloud cover (fraction)
                case(17)
                output(:, 1) = cloudc
                l_3d = .false.

                ! clstr  = stratiform cloud cover (fraction)
                case(18)
                output(:, 1) = clstr
                l_3d = .false.

                ! cltop  = norm. pressure at cloud top
                case(19)
                output(:, 1) = cltop
                l_3d = .false.

                ! prtop  = top of precipitation (level index)
                case(20)
                output(:, 1) = prtop
                l_3d = .false.

                ! precnv = convective precipitation  [g/(m^2 s)], total
                case(31)
                output(:, 1) = precnv
                l_3d = .false.

                ! precls = large-scale precipitation [g/(m^2 s)], total
                case(32)
                output(:, 1) = precls / 3600.0_dp
                l_3d = .false.

                ! snowcv = convective precipitation  [g/(m^2 s)], snow only
                !case(33)
                !output(:, 1) = snowcv

                ! snowls = large-scale precipitation [g/(m^2 s)], snow only
                !case(34)
                !output(:, 1) = snowls

                ! cbmf   = cloud-base mass flux
                case(35)
                output(:, 1) = cbmf
                l_3d = .false.

                ! tsr    = top-of-atm. shortwave radiation (downward)
                case(36)
                output(:, 1) = tsr
                l_3d = .false.

                ! ssrd   = surface shortwave radiation (downward-only)
                case(37)
                output(:, 1) = ssrd
                l_3d = .false.

                ! ssr    = surface shortwave radiation (net downward)
                case(38)
                output(:, 1) = ssr
                l_3d = .false.

                ! slrd   = surface longwave radiation  (downward-only)
                case(39)
                output(:, 1) = slrd
                l_3d = .false.

                ! slr    = surface longwave radiation  (net upward)
                case(40)
                output(:, 1) = slr
                l_3d = .false.

                ! olr    = outgoing longwave radiation (upward)
                case(41)
                output(:, 1) = olr
                l_3d = .false.

                ! slru   = surface longwave emission   (upward)
                !                                   (1:land, 2:sea, 3: wgt. average)
                case(42)
                output(:, 1) = slru(:, 3)
                l_3d = .false.

                ! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
                case(43)
                output(:, 1) = ustr(:, 3)
                l_3d = .false.

                ! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
                case(44)
                output(:, 1) = vstr(:, 3)
                l_3d = .false.

                ! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
                case(45)
                output(:, 1) = shf(:, 3)
                l_3d = .false.

                ! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
                case(46)
                output(:, 1) = evap(:, 3)
                l_3d = .false.

                ! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
                case(47)
                output(:, 1) = hfluxn(:, 3)
                l_3d = .false.

                ! alb_l = surface albedo over land
                case(48)
                output(:, 1) = alb_l
                l_3d = .false.

                ! alb_s = surface albedo over sea
                case(49)
                output(:, 1) = alb_s
                l_3d = .false.

                ! albsfc = surface albedo
                case(50)
                output(:, 1) = albsfc
                l_3d = .false.

                ! snowc = snow cover
                case(51)
                output(:, 1) = snowc
                l_3d = .false.

                ! fsol = flux of incoming solar radiation
                case(52)
                output(:, 1) = fsol
                l_3d = .false.

                ! ozone = ozone concentration
                case(53)
                output(:, 1) = ozone
                l_3d = .false.

                ! ozupp = ozone concentration ...
                case(54)
                output(:, 1) = ozupp
                l_3d = .false.

                ! zenit = ???
                case(55)
                output(:, 1) = zenit
                l_3d = .false.

                ! stratz = ???
                case(56)
                output(:, 1) = stratz
                l_3d = .false.

                ! stl_am = Surface temperature from the land model
                case(57)
                output(:, 1) = stl_am
                l_3d = .false.

                ! soilw_am = Soil water from the land model
                case(58)
                output(:, 1) = soilw_am
                l_3d = .false.

                ! sst_am = Sea surface temperature from the ocean model
                case(59)
                output(:, 1) = sst_am
                l_3d = .false.

                ! ssti_om = Sea surface + ice temperature from the ocean model
                case(60)
                output(:, 1) = ssti_om
                l_3d = .false.

                ! phis0 = spectrally filtered surface geopotential
                case(61)
                output(:, 1) = pack(phis0, .true.)
                l_3d = .false.

                ! fmask1 = fractional land-sea mask
                case(62)
                output(:, 1) = pack(fmask1, .true.)
                l_3d = .false.

                ! tt_cnv  =  temperature tendency due to convection
                case(101)
                output = tt_cnv

                ! qt_cnv  = sp. humidity tendency due to convection
                case(102)
                output = qt_cnv

                ! tt_lsc  =  temperature tendency due to large-scale condensation
                case(103)
                output = tt_lsc

                ! qt_lsc  = sp. humidity tendency due to large-scale condensation
                case(104)
                output = qt_lsc

                ! tt_rsw  =  temperature tendency due to short-wave radiation
                case(105)
                output = tt_rsw

                ! tt_rlw  =  temperature tendency due to long-wave radiation
                case(106)
                output = tt_rlw

                ! ut_sflx  =       u-wind tendency due to surface fluxes
                case(107)
                output = ut_sflx

                ! vt_sflx  =       v-wind tendency due to surface fluxes
                case(108)
                output = vt_sflx

                ! tt_sflx  =  temperature tendency due to surface fluxes
                case(109)
                output = tt_sflx

                ! qt_sflx  = sp. humidity tendency due to surface fluxes
                case(110)
                output = qt_sflx

                ! ut_pbl  =       u-wind tendency due to PBL and diffusive processes
                case(111)
                output = ut_pbl

                ! vt_pbl  =       v-wind tendency due to PBL and diffusive processes
                case(112)
                output = vt_pbl

                ! tt_pbl  =  temperature tendency due to PBL and diffusive processes
                case(113)
                output = tt_pbl

                ! qt_pbl  = sp. humidity tendency due to PBL and diffusive processes
                case(114)
                output = qt_pbl

                ! ut_phy  =       u-wind tendency due to all physics processes
                case(115)
                output = ut_phy

                ! vt_phy  =       v-wind tendency due to all physics processes
                case(116)
                output = vt_phy

                ! tt_phy  =  temperature tendency due to all physics processes
                case(117)
                output = tt_phy

                ! qt_phy  = sp. humidity tendency due to all physics processes
                case(118)
                output = qt_phy

                ! ut_sppt =       u-wind tendency due to stochastic perturbation
                case(119)
                output = ut_sppt

                ! vt_sppt =       v-wind tendency due to stochastic perturbation
                case(120)
                output = vt_sppt

                ! tt_sppt =  temperature tendency due to stochastic perturbation
                case(121)
                output = tt_sppt

                ! qt_sppt = sp. humidity tendency due to stochastic perturbation
                case(122)
                output = qt_sppt

                ! 3D Stochastic perturbation pattern
                case(123)
                output = sppt

                case default
                print *, 'Variable no.', varID, ' unavailable for output'
            end select
        end subroutine fetch_output_grid

        subroutine add_var_info(varID, file_ID, nc_var_ID, dimids)
            integer, intent(in) :: varID, file_ID
            integer, intent(inout) :: nc_var_ID
            integer, intent(in) :: dimids(4)
            integer :: dimids_2d(3)
            logical :: l_3d
            character(len=128) :: name
            character(len=32) :: units

            l_3d = .true.

            select case(varID)
                ! ug1/vg1 (ms-1)
                case(1)
                name = 'zonal_velocity'
                units = 'm s-1'

                 ! vg1
                case(2)
                name = 'meridional_velocity'
                units = 'm s-1'

                ! tg1 (K)
                case(3)
                name = 'temperature'
                units = 'K'

                ! qg1 (kg/kg)
                case(4)
                name = 'specific_humidity'
                units = ''

                ! phig1  = geopotential
                case(5)
                name = 'geopotential_height'
                units = 'm'

                ! pslg1
                case(6)
                name = 'logarithm_of_surface_pressure'
                units = ''
                l_3d = .false.

                ! se     = dry static energy
                case(7)
                name = 'dry_static_energy'
                units = 'J kg-1'

                ! rh     = relative humidity
                case(8)
                name = 'relative_humidity'
                units = ''

                ! qsat   = saturation specific humidity (g/kg)
                case(9)
                name = 'saturation_specific_humidity'
                units = 'g kg-1'

                ! psg    = surface pressure
                case(10)
                name = 'surface_pressure'
                units = 'Pa'
                l_3d = .false.

                ! ts     = surface temperature
                case(11)
                name = 'surface_temperature'
                units = 'K'
                l_3d = .false.

                ! tskin  = skin temperature
                case(12)
                name = 'skin_temperature'
                units = 'K'
                l_3d = .false.

                ! u0     = near-surface u-wind
                case(13)
                name = 'near_surface_zonal_velocity'
                units = 'm s-1'
                l_3d = .false.

                ! v0     = near-surface v-wind
                case(14)
                name = 'near_surface_meridional_velocity'
                units = 'm s-1'
                l_3d = .false.

                ! t0     = near-surface air temperature
                case(15)
                name = 'near_surface_temperature'
                units = 'K'
                l_3d = .false.

                ! q0     = near-surface specific humidity (g/kg)
                case(16)
                name = 'near_surface_specific_humidity'
                units = 'g kg-1'
                l_3d = .false.

                ! cloudc = total cloud cover (fraction)
                case(17)
                name = 'cloud_cover'
                units = ''
                l_3d = .false.

                ! clstr  = stratiform cloud cover (fraction)
                case(18)
                name = 'stratiform_cloud_cover'
                units = ''
                l_3d = .false.

                ! cltop  = norm. pressure at cloud top
                case(19)
                name = 'cloud_top_pressure'
                units = ''
                l_3d = .false.

                ! prtop  = top of precipitation (level index)
                case(20)
                name = 'level_of_precipitation'
                units = ''
                l_3d = .false.

                ! precnv = convective precipitation  [g/(m^2 s)], total
                case(31)
                name = 'convective_precipitation'
                units = 'g m-2 s-1'
                l_3d = .false.

                ! precls = large-scale precipitation [g/(m^2 s)], total
                case(32)
                name = 'large_scale_precipitation'
                units = 'g m-2 s-1'
                l_3d = .false.

                ! snowcv = convective precipitation  [g/(m^2 s)], snow only
                !case(33)

                ! snowls = large-scale precipitation [g/(m^2 s)], snow only
                !case(34)

                ! cbmf   = cloud-base mass flux
                case(35)
                name = 'cloud_base_mass_flux'
                units ='unknown'
                l_3d = .false.

                ! tsr    = top-of-atm. shortwave radiation (downward)
                case(36)
                name = 'top_of_atmosphere_shortwave_radiation'
                units = 'unknown'
                l_3d = .false.

                ! ssrd   = surface shortwave radiation (downward-only)
                case(37)
                name = 'downward_shortwave_radiation_at_surface'
                units = 'unknown'
                l_3d = .false.

                ! ssr    = surface shortwave radiation (net downward)
                case(38)
                name = 'net_downward_shortwave_radiation_at_surface'
                units = 'unknown'
                l_3d = .false.

                ! slrd   = surface longwave radiation  (downward-only)
                case(39)
                name = 'downward_longwave_radiation_at_surface'
                units = 'unknown'
                l_3d = .false.

                ! slr    = surface longwave radiation  (net upward)
                case(40)
                name = 'net_upward_longwave_radiation_at_surface'
                units = 'unknown'
                l_3d = .false.

                ! olr    = outgoing longwave radiation (upward)
                case(41)
                name = 'outgoing_longwave_radiation'
                units = 'unknown'
                l_3d = .false.

                ! slru   = surface longwave emission   (upward)
                !                                   (1:land, 2:sea, 3: wgt. average)
                case(42)
                name = 'surface_longwave_emission'
                units = 'unknown'
                l_3d = .false.

                ! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
                case(43)
                name = 'zonal_wind_stress'
                units = 'unknown'
                l_3d = .false.

                ! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
                case(44)
                name = 'meridional_wind_stress'
                units = 'unknown'
                l_3d = .false.

                ! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
                case(45)
                name = 'sensible_heat_flux'
                units = 'unknown'
                l_3d = .false.

                ! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
                case(46)
                name = 'evaporation'
                units = 'g m-2 s-1'
                l_3d = .false.

                ! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
                case(47)
                name = 'net_heat_flux_into_surface'
                units = ''
                l_3d = .false.

                ! alb_l = Surface albedo over land
                case(48)
                name = 'surface_albedo_over_land'
                units = 'unknown'
                l_3d = .false.

                ! alb_s = surface albedo over sea
                case(49)
                name = 'surface_albedo_over_sea'
                units = 'unknown'
                l_3d = .false.

                ! albsfc = surface albedo
                case(50)
                name = 'surface_albedo'
                units = 'unknown'
                l_3d = .false.

                ! snowc = snow cover
                case(51)
                name = 'snow_cover'
                units = 'unknown'
                l_3d = .false.

                ! fsol = flux of incoming solar radiation
                case(52)
                name = 'flux_of_incoming_solar_radiation'
                units = 'unknown'
                l_3d = .false.

                ! ozone = ozone concentration
                case(53)
                name = 'ozone'
                units = 'unknown'
                l_3d = .false.

                ! ozupp = ozone concentration ...
                case(54)
                name = 'ozupp'
                units = 'unknown'
                l_3d = .false.

                ! zenit = ???
                case(55)
                name = 'zenit'
                units = 'unknown'
                l_3d = .false.

                ! stratz = ???
                case(56)
                name = 'stratz'
                units = 'unknown'
                l_3d = .false.

                ! stl_am = Surface temperature from the land model
                case(57)
                name = 'land_model_surface_temperature'
                units = 'K'
                l_3d = .false.

                ! soilw_am = Soil water from the land model
                case(58)
                name = 'land_model_soil_moisture'
                units = 'unknown'
                l_3d = .false.

                ! sst_am = Sea surface temperature from the ocean model
                case(59)
                name = 'ocean_model_surface_temperature'
                units = 'unknown'
                l_3d = .false.

                ! ssti_om = Sea surface + ice temperature from the ocean model
                case(60)
                name = 'ocean_and_ice_model_surface_temperature'
                units = 'K'
                l_3d = .false.

                ! phis0 = spectrally filtered surface geopotential
                case(61)
                name = 'surface_geopotential'
                units = 'unknown'
                l_3d = .false.

                ! fmask1 = fractional land-sea mask
                case(62)
                name = 'fractional_land_sea_mask'
                units = ''
                l_3d = .false.

                ! tt_cnv  =  temperature tendency due to convection
                case(101)
                name = 'temperature_tendency_due_to_convection'
                units = 'K s-1'

                ! qt_cnv  = sp. humidity tendency due to convection
                case(102)
                name = 'specific_humidity_tendency_due_to_convection'
                units = 's-1'

                ! tt_lsc  =  temperature tendency due to large-scale condensation
                case(103)
                name = 'temperature_tendency_due_to_condensation'
                units = 'K s-1'

                ! qt_lsc  = sp. humidity tendency due to large-scale condensation
                case(104)
                name = 'specific_humidity_tendency_due_to_condensation'
                units = 's-1'

                ! tt_rsw  =  temperature tendency due to short-wave radiation
                case(105)
                name = 'temperature_tendency_due_to_shortwave_radiation'
                units = 'K s-1'

                ! tt_rlw  =  temperature tendency due to long-wave radiation
                case(106)
                name = 'temperature_tendency_due_to_longwave_radiation'
                units = 'K s-1'

                ! ut_sflx  =       u-wind tendency due to surface fluxes
                case(107)
                name = 'zonal_velocity_tendency_due_to_surface_fluxes'
                units = 'm s-2'

                ! vt_sflx  =       v-wind tendency due to surface fluxes
                case(108)
                name = 'meridional_velocity_tendency_due_to_surface_fluxes'
                units = 'm s-2'

                ! tt_sflx  =  temperature tendency due to surface fluxes
                case(109)
                name = 'temperature_tendency_due_to_surface_fluxes'
                units = 'K s-1'

                ! qt_sflx  = sp. humidity tendency due to surface fluxes
                case(110)
                name = 'specific_humidity_tendency_due_to_surface_fluxes'
                units = 's-1'

                ! ut_pbl  =       u-wind tendency due to PBL and diffusive processes
                case(111)
                name = 'zonal_velocity_tendency_due_to_vertical_diffusion'
                units = 'm s-2'

                ! vt_pbl  =       v-wind tendency due to PBL and diffusive processes
                case(112)
                name = 'meridional_velocity_tendency_due_to_vertical_diffusion'
                units = 'm s-2'

                ! tt_pbl  =  temperature tendency due to PBL and diffusive processes
                case(113)
                name = 'temperature_tendency_due_to_vertical_diffusion'
                units = 'K s-1'

                ! qt_pbl  = sp. humidity tendency due to PBL and diffusive processes
                case(114)
                name = 'specific_humidity_tendency_due_to_vertical_diffusion'
                units = 's-1'

                ! ut_phy  =       u-wind tendency due to all physics processes
                case(115)
                name = 'zonal_velocity_tendency_due_to_all_parametrizations'
                units = 'm s-2'

                ! vt_phy  =       v-wind tendency due to all physics processes
                case(116)
                name = 'meridional_velocity_tendency_due_to_all_parametrizations'
                units = 'm s-2'

                ! tt_phy  =  temperature tendency due to all physics processes
                case(117)
                name = 'temperature_tendency_due_to_all_parametrizations'
                units = 'K s-1'

                ! qt_phy  = sp. humidity tendency due to all physics processes
                case(118)
                name = 'specific_humidity_tendency_due_to_all_parametrizations'
                units = 's-1'

                ! ut_sppt =       u-wind tendency due to stochastic perturbation
                case(119)
                name = 'zonal_velocity_tendency_due_to_stochastic_perturbation'
                units = 'm s-2'

                ! vt_sppt =       v-wind tendency due to stochastic perturbation
                case(120)
                name = 'meridional_velocity_tendency_due_to_stochastic_perturbation'
                units = 'm s-2'

                ! tt_sppt =  temperature tendency due to stochastic perturbation
                case(121)
                name = 'temperature_tendency_due_to_stochastic_perturbation'
                units = 'K s-1'

                ! qt_sppt = sp. humidity tendency due to stochastic perturbation
                case(122)
                name = 'specific_humidity_tendency_due_to_stochastic_perturbation'
                units = 's-1'

                ! 3D Stochastic perturbation pattern
                case(123)
                name = 'stochastic_perturbation'
                units = ''

                case default
                print *, 'Variable no.', varID, ' unavailable for output'
            end select

            ! Don't include the vertical dimension on 2d fields
            if (l_3d) then
                call check( nf90_def_var(file_ID, trim(name), NF90_DOUBLE, dimids, nc_var_ID) )
            else
                dimids_2d(1) = dimids(1)
                dimids_2d(2) = dimids(2)
                dimids_2d(3) = dimids(4)
                call check( nf90_def_var(file_ID, trim(name), NF90_DOUBLE, dimids_2d, nc_var_ID) )
            end if

            call check( nf90_put_att(file_ID, nc_var_ID, 'units', trim(units)) )
        end subroutine add_var_info

        subroutine init_nc(stream)
            type(output_stream), intent(inout) :: stream
            integer :: lon_dimid, lat_dimid, lvl_dimid, rec_dimid, &
                    lon_varid, lat_varid, lvl_varid
            integer :: dimids(4)
            integer :: n

            ! Open the netcdf file, assigning the file_ID
            call check( nf90_create(trim(stream%filename), NF90_CLOBBER, stream%file_ID) )

            ! Add the dimensions of the model to the netcdf file
            call check( nf90_def_dim(stream%file_ID, 'longitude', ix, lon_dimid) )
            call check( nf90_def_dim(stream%file_ID, 'latitude' , il, lat_dimid) )
            if (stream%plevs) then
                call check( nf90_def_dim(stream%file_ID, 'pressure', kx, lvl_dimid) )
            else
                call check( nf90_def_dim(stream%file_ID, 'sigma'   , kx, lvl_dimid) )
            end if
            call check( nf90_def_dim(stream%file_ID, 'forecast_period', NF90_UNLIMITED, rec_dimid) )

            ! Add the dimensions as variables as well
            call check( nf90_def_var(stream%file_ID, 'longitude', NF90_DOUBLE, lon_dimid, lon_varid) )
            call check( nf90_def_var(stream%file_ID, 'latitude' , NF90_DOUBLE, lat_dimid, lat_varid) )
            if (stream%plevs) then
                call check( nf90_def_var(stream%file_ID, 'pressure' , NF90_DOUBLE, lvl_dimid, lvl_varid) )
            else
                call check( nf90_def_var(stream%file_ID, 'sigma'    , NF90_DOUBLE, lvl_dimid, lvl_varid) )
            end if
            call check( nf90_def_var(stream%file_ID, 'forecast_period' , NF90_DOUBLE, rec_dimid, stream%rec_varid) )

            ! Assign units attributes to coordinate variables.
            call check( nf90_put_att(stream%file_ID, lon_varid, 'units', 'degrees') )
            call check( nf90_put_att(stream%file_ID, lat_varid, 'units', 'degrees') )
            if (stream%plevs) then
                call check( nf90_put_att(stream%file_ID, lvl_varid, 'units', 'hPa') )
            end if
            call check( nf90_put_att(stream%file_ID, stream%rec_varid, 'units', 'seconds') )

            ! Define the netCDF variables for the output fields
            dimids = (/ lon_dimid, lat_dimid, lvl_dimid, rec_dimid /)
            do n=1, stream%nvars
                call add_var_info(stream%var_ID(n), stream%file_ID, stream%nc_var_ID(n), dimids)
            end do

            ! End define mode.
            call check( nf90_enddef(stream%file_ID) )

            ! Write the coordinate variable data
            call check( nf90_put_var(stream%file_ID, lon_varid, (/ (n*(360.0_dp/ix), n=0, ix-1) /)) )
            call check( nf90_put_var(stream%file_ID, lat_varid, deglat_s(il:1:-1)%val) )
            if (stream%plevs) then
                call check( nf90_put_var(stream%file_ID, lvl_varid, pout(kx:1:-1)*1000) )
            else
                call check( nf90_put_var(stream%file_ID, lvl_varid, sig(kx:1:-1)%val) )
            end if
            call check( nf90_put_var(stream%file_ID, stream%rec_varid, 0.0_dp) )
        end subroutine init_nc

        subroutine check(status)
            ! Wrapper subroutine for NetCDF library functions that checks the returned
            ! status and stops the program on any errors

            ! Status identifier output from netCDF library functions
            integer, intent (in) :: status

            if(status /= NF90_NOERR) then
                print *, NF90_STRERROR(status)
                stop "Stopped"
            end if
        end subroutine check
end module
