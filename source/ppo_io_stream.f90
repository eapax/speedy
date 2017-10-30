module ppo_IO_stream

    use mod_atparam
    use mod_dynvar
    use mod_physvar

    implicit none

    private
    public IO_stream, init_IO_stream, write_IO_stream

    integer, parameter :: recl_spec = 4*mx*nx
    integer, parameter :: recl_grid = 4*ix*il

    ! Generic type used to describe an output stream
    type IO_stream
        ! Spectral and Grid variables need to be treated differently in
        ! separate streams due to their different grids and the use of complex
        ! numbers
        logical :: spectral

        ! The ID used in open/read/write statements
        integer :: file_ID

        ! An array of integers corresponding to the variables being output
        integer, allocatable :: var_ID(:)

        ! The record counter used for direct access in read/write statements
        integer :: rec = 1
    end type IO_stream

    contains
        !
        function init_IO_stream(filename, spectral, var_ID) result(stream)
            character(len=*), intent(in) :: filename
            logical, intent(in) :: spectral
            integer, intent(in) :: var_ID(:)
            type(IO_stream) :: stream

            integer :: recl

            stream%spectral = spectral
            stream%file_ID = 200
            allocate(stream%var_ID(size(var_ID)))
            stream%var_ID = var_ID

            if (spectral) then
                recl = recl_spec
            else
                recl = recl_grid
            end if

            open(stream%file_ID, file=filename, form='unformatted', &
                    access='direct', recl=recl)
        end function

        ! Write each of the variables associated with the given stream to file.
        ! Acts a wrapper to separate write functions for spectral and grid
        ! variables.
        subroutine write_IO_stream(stream)
            type(IO_stream), intent(inout) :: stream

            if (stream%spectral) then
                call write_spectral(stream)
            else
                call write_grid(stream)
            end if
        end subroutine

        subroutine write_spectral(stream)
            type(IO_stream), intent(inout) :: stream
            complex :: output(mx, nx, kx)
            real(4) :: re_output(mx, nx, kx), im_output(mx, nx, kx)
            integer :: n, k

            do n=1, size(stream%var_ID)
                output = fetch_output_spectral(stream%var_ID(n))

                ! Write real and imaginary parts as separate variables
                re_output = real(real(output))
                im_output = real(aimag(output))

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
            type(IO_stream), intent(inout) :: stream
            real(4) :: output(ix*il, kx)
            integer :: n, k

            do n=1, size(stream%var_ID)
                output = fetch_output(stream%var_ID(n))

                ! For some reason the height levels need to be written backwards
                do k=kx, 1, -1
                    write(stream%file_ID, rec=stream%rec) output(:, k)
                    stream%rec = stream%rec + 1
                end do
            end do
        end subroutine write_grid

        ! Get the variable corresponding to the varID in single bit precision
        ! Essentially a look up for all the variables in mod_dynvar
        ! TODO implement 2d variables in same interface
        function fetch_output_spectral(varID) result(output)
            integer :: varID
            complex :: output(mx, nx, kx)

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
        function fetch_output(varID) result(output)
            integer :: varID
            real(4) :: output(ix*il, kx)

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
                output = qg1

                ! phig1  = geopotential
                case(5)
                output=phig1

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
                case(21)
                output = tt_cnv

                ! qt_cnv  = sp. humidity tendency due to convection
                case(22)
                output = qt_cnv

                ! tt_lsc  =  temperature tendency due to large-scale condensation
                case(23)
                output = tt_lsc

                ! qt_lsc  = sp. humidity tendency due to large-scale condensation
                case(24)
                output = qt_lsc

                ! tt_rsw  =  temperature tendency due to short-wave radiation
                case(25)
                output = tt_rsw

                ! tt_rlw  =  temperature tendency due to long-wave radiation
                case(26)
                output = tt_rlw

                ! ut_pbl  =       u-wind tendency due to PBL and diffusive processes
                case(27)
                output = ut_pbl

                ! vt_pbl  =       v-wind tendency due to PBL and diffusive processes
                case(28)
                output = vt_pbl

                ! tt_pbl  =  temperature tendency due to PBL and diffusive processes
                case(29)
                output = tt_pbl

                ! qt_pbl  = sp. humidity tendency due to PBL and diffusive processes
                case(30)
                output = qt_pbl

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