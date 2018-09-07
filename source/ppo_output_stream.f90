module ppo_output_stream

    use mod_atparam
    use mod_dynvar
    use mod_physvar
    use mod_physcon, only: gg
    use mod_date, only: imonth, month_start
    use spectral, only: uvspec, grid
    use ppo_plevs, only: pressure_levels, np
    use mod_prec, only: sp, dp

    implicit none

    private
    public initialise_output, update_output

    integer :: recl_spec, recl_grid

    ! Counter so that each output stream has a unique file ID
    integer :: next_file_ID = 200

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

        ! An array of integers corresponding to the variables being output
        integer :: nvars
        integer, allocatable :: var_ID(:)

        ! The record counter used for direct access in read/write statements
        integer :: recl
        integer :: rec = 1
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

            ! Set record lengths to grid size
            recl_spec = 4*mx*nx
            recl_grid = 4*ngp

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
            end if

            ! Temperature
            if (l_transform_field(3) .or. l_transform_field(5)) then
                do k=1,kx
                    call grid(t(:,:,k,1), tg1(:,k), 1)
                end do
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

        !
        function init_output_stream(filename, spectral, plevs, nstpinc, &
                nstpout, nstpopen, nvars, var_ID) &
                result(stream)
            character(len=*), intent(in) :: filename
            logical, intent(in) :: spectral, plevs
            integer, intent(in) :: nstpinc, nstpout, nstpopen
            integer, intent(in) :: nvars, var_ID(:)
            type(output_stream) :: stream

            integer :: recl

            ! Pass parameters to the new object
            stream%filename = filename
            stream%spectral = spectral
            stream%plevs = plevs
            stream%nstpinc = nstpinc
            stream%nstpout = nstpout
            stream%nstpopen = nstpopen
            stream%nvars = nvars
            allocate(stream%var_ID(stream%nvars))
            stream%var_ID = var_ID

            if (spectral) then
                stream%recl = recl_spec
            else
                stream%recl = recl_grid
            end if

            ! Assign a unique file ID to the output stream
            stream%file_ID = next_file_ID
            next_file_ID = next_file_ID + 1
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

            if (xmod(istep, stream%nstpopen)) call reinit_output_stream(stream)

            if (xmod(istep, stream%nstpinc)) call incr_output_stream(stream)

            if (xmod(istep, stream%nstpout)) call write_output_stream(stream)
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

        ! todo Close the currently open file and open a new file named by the
        ! current time
        subroutine reinit_output_stream(stream)
            type(output_stream), intent(inout) :: stream

            close(stream%file_ID)

            open(stream%file_ID, file=trim(stream%filename), &
                    form='unformatted', access='direct', recl=stream%recl)
        end subroutine

        ! todo Increment the output on substeps relative to how often it is
        ! written to file
        subroutine incr_output_stream(stream)
            type(output_stream), intent(inout) :: stream
        end subroutine

        ! Write each of the variables associated with the given stream to file.
        ! Acts a wrapper to separate write functions for spectral and grid
        ! variables.
        subroutine write_output_stream(stream)
            type(output_stream), intent(inout) :: stream

            if (stream%spectral) then
                call write_spectral(stream)
            else
                call write_grid(stream)
            end if
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
            real(sp) :: output_sp(ngp, kx)
            real(sp) :: output_p(ngp, np)
            integer :: n, k

            do n=1, stream%nvars
                output = fetch_output_grid(stream%var_ID(n))

                if(stream%plevs) then
                    ! Interpolate output to pressure levels
                    call pressure_levels(stream%var_ID(n), output, output_p)
                    do k=np, 1, -1
                        write(stream%file_ID, rec=stream%rec) output_p(:, k)
                        stream%rec = stream%rec + 1
                    end do
                else
                    ! Otherwise write model level output
                    if (stream%var_ID(n)==5) then
                        output_sp = output / gg !m
                    else
                        output_sp = output
                    end if

                    ! For some reason the height levels need to be written backwards
                    do k=kx, 1, -1
                        write(stream%file_ID, rec=stream%rec) output_sp(:, k)
                        stream%rec = stream%rec + 1
                    end do
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
        ! TODO implement 2d variables in same interface
        function fetch_output_grid(varID) result(output)
            integer :: varID
            real(dp) :: output(ngp, kx)

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
                !case(6)
                !output(:, 1) = pslg1

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
                !case(10)
                !output(:, 1) = psg

                ! ts     = surface temperature
                !case(11)
                !output(:, 1) = ts

                ! tskin  = skin temperature
                !case(12)
                !output(:, 1) = tskin

                ! u0     = near-surface u-wind
                !case(13)
                !output(:, 1) = u0

                ! v0     = near-surface v-wind
                !case(14)
                !output(:, 1) = v0

                ! t0     = near-surface air temperature
                !case(15)
                !output(:, 1) = t0

                ! q0     = near-surface specific humidity (g/kg)
                !case(16)
                !output(:, 1) = q0

                ! cloudc = total cloud cover (fraction)
                !case(17)
                !output(:, 1) = cloudc

                ! clstr  = stratiform cloud cover (fraction)
                !case(18)
                !output(:, 1) = clstr

                ! cltop  = norm. pressure at cloud top
                !case(19)
                !output(:, 1) = cltop

                ! prtop  = top of precipitation (level index)
                !case(20)
                !output(:, 1) = prtop

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

                ! precnv = convective precipitation  [g/(m^2 s)], total
                !case(31)
                !output(:, 1) = precnv

                ! precls = large-scale precipitation [g/(m^2 s)], total
                !case(32)
                !output(:, 1) = precls

                ! snowcv = convective precipitation  [g/(m^2 s)], snow only
                !case(33)
                !output(:, 1) = snowcv

                ! snowls = large-scale precipitation [g/(m^2 s)], snow only
                !case(34)
                !output(:, 1) = snowls

                ! cbmf   = cloud-base mass flux
                !case(35)
                !output(:, 1) = cbmf

                ! tsr    = top-of-atm. shortwave radiation (downward)
                !case(36)
                !output(:, 1) = tsr

                ! ssrd   = surface shortwave radiation (downward-only)
                !case(37)
                !output(:, 1) = ssrd

                ! ssr    = surface shortwave radiation (net downward)
                !case(37)
                !output(:, 1) = ssr

                ! slrd   = surface longwave radiation  (downward-only)
                !case(37)
                !output(:, 1) = slrd

                ! slr    = surface longwave radiation  (net upward)
                !case(37)
                !output(:, 1) = slr

                ! olr    = outgoing longwave radiation (upward)
                !case(37)
                !output(:, 1) = olr

                ! slru   = surface longwave emission   (upward)
                !                                   (1:land, 2:sea, 3: wgt. average)
                !case(37)
                !output(:, 1:3) = slru

                ! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
                !case(37)
                !output(:, 1:3) = ustr

                ! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
                !case(37)
                !output(:, 1:3) = vstr

                ! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
                !case(37)
                !output(:, 1:3) = shf

                ! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
                !case(37)
                !output(:, 1:3) = evap

                ! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
                !case(37)
                !output(:, 1:3) = hfluxn

                case default
                print *, 'Variable no.', varID, ' unavailable for output'
            end select
        end function
end module