!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> A module for computing SPPT patterns to be used as multiplicative noise applied to physical tendencies
!> Stochastically Perturbed Parametrization Tendencies (SPPT) is a parametrization of model error.
!> See ECMWF Tech. Memo. #598 (Palmer et al. 2009)
module phy_sppt
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public sppt_on, setup_sppt, ini_sppt, truncate_sppt, gen_sppt

    namelist /sppt/ sppt_on, nscales
    namelist /sppt_parameters/ mu, time_decorr, len_decorr, stddev

    ! Turn on SPPT?
    logical :: sppt_on = .false.

    ! Array for tapering value of SPPT in the different layers of the atmosphere
    ! A value of 1 means the tendency is not tapered at that level
    type(rpe_var), allocatable :: mu(:)

    ! Number of correlation scales for SPPT perturbations
    integer :: nscales = 3

    ! Namelist parameters used to setup SPPT constants
    ! Decorrelation time of SPPT perturbation (in hours)
    real(dp), allocatable :: time_decorr(:)
    ! Correlation length scale of SPPT perturbation (in metres)
    real(dp), allocatable :: len_decorr(:)
    ! Standard deviation of SPPT perturbation (in grid point space)
    real(dp), allocatable :: stddev(:)

    ! SPPT parameters initialised in ini_sppt
    ! Time autocorrelation of spectral AR(1) signals
    type(rpe_var), allocatable :: phi(:)
    ! Total wavenumber-wise standard deviation of spectral signals
    type(rpe_var), allocatable :: sigma(:, :, :)

    ! Perturbations in spectral space
    type(rpe_complex_var), allocatable :: sppt_spec(:, :, :)

    logical :: first = .true.

    contains
        subroutine setup_sppt(fid)
            integer, intent(in) :: fid

            read(fid, sppt)
            write(*, sppt)

            allocate(mu(kx))
            allocate(time_decorr(nscales))
            allocate(len_decorr(nscales))
            allocate(stddev(nscales))
            allocate(phi(nscales))
            allocate(sigma(mx, nx, nscales))
            allocate(sppt_spec(mx, nx, nscales))

            read(fid, sppt_parameters)
            write(*, sppt_parameters)
        end subroutine setup_sppt

        subroutine truncate_sppt()
            call apply_truncation(mu)
            call apply_truncation(phi)
            call apply_truncation(sigma)
        end subroutine truncate_sppt

        !> @brief
        !> Generate grid point space SPPT pattern
        !> distribution.
        !> @return sppt_grid the generated grid point pattern
        function gen_sppt() result(sppt_grid_out)
            use spectral, only: grid

            integer :: m, n, k
            type(rpe_var) :: sppt_grid(ix, il), sppt_grid_total(ix, il), &
                             sppt_grid_out(ngp, kx)
            type(rpe_complex_var) :: eta(mx, nx, nscales)
            type(rpe_var) :: randreal, randimag

            ! Generate Gaussian noise
            do k=1, nscales
                do n=1, nx
                    do m=1, mx
                        randreal = randn( &
                                rpe_literal(0.0_dp), rpe_literal(1.0_dp))
                        randimag = randn( &
                                rpe_literal(0.0_dp), rpe_literal(1.0_dp))

                        ! Clip noise to +- 10 standard deviations
                        eta(m,n,k) = cmplx(&
                                min(rpe_literal(10.0_dp), abs(randreal)) * &
                                        sign(rpe_literal(1.0_dp), randreal),&
                                min(rpe_literal(10.0_dp), abs(randimag)) * &
                                        sign(rpe_literal(1.0_dp), randimag))
                    end do
                end do
            end do

            ! If first timestep
            if (first) then
                ! First AR(1) step
                do n=1, nscales
                    sppt_spec(:,:,n) = &
                            (1 - phi(n)**2)**(rpe_literal(-0.5_dp)) * &
                                    sigma(:,:,n) * eta(:,:,n)
                end do
                first = .false.
            else
                ! Subsequent AR(1) steps
                do n=1, nscales
                    sppt_spec(:,:,n) = phi(n)*sppt_spec(:,:,n) + &
                            sigma(:,:,n) * eta(:,:,n)
                end do
            end if

            ! Sum SPPT perturbations over correlation scales
            sppt_grid_total = 0.0_dp
            do n=1, nscales
                ! Convert to grid point space
                call grid(sppt_spec(:,:,n), sppt_grid, 1)
                sppt_grid_total = sppt_grid_total + sppt_grid
            end do

            ! Clip to +/- 1.0
            sppt_grid_total = min(rpe_literal(1.0_dp), abs(sppt_grid_total)) * &
                    sign(rpe_literal(1.0_dp), sppt_grid_total)

            ! Apply tapering
            do k=1,kx
                sppt_grid_out(:, k) = rpe_literal(1.0_dp) + &
                        reshape(sppt_grid_total, (/ngp/)) * mu(k)
            end do
        end function gen_sppt

        !> @brief
        !> Calculates the phi and sigma parameters in the SPPT calculations
        !> from the defined correlation scales and standard deviations
        subroutine ini_sppt()
            use mod_tsteps, only: nsteps
            use mod_dyncon1, only: rearth
            use spectral, only: el2

            integer :: n, sc
            real(dp) :: f0(nscales)

            ! Calculate time autocorrelation factor as a function of timestep
            phi = exp(-(24.0_dp/real(nsteps, dp))/time_decorr)

            ! Generate spatial amplitude pattern
            do sc=1, nscales
                f0(sc) = sum((/ ((2*n+1)*exp( &
                        -0.5_dp*(len_decorr(sc)/rearth)**2*n*(n+1) &
                        ), n=1,ntrun) /))
                f0(sc) = sqrt((stddev(sc)**2*(1-phi(sc)**2))/f0(sc))
                sigma(:,:,sc) = f0(sc) * exp(-0.25_dp*len_decorr(sc)**2 * el2)
            end do

            ! Initialise the random number generator
            call time_seed()
        end subroutine ini_sppt

        !> @brief
        !> Generates a random number drawn for the specified normal
        !> distribution.
        !> @param mean the mean of the distribution to draw from
        !> @param stdev the standard deviation of the distribution to draw from
        !> @return randn the generated random number
        function randn(mean, stdev)
            use mod_prec

            type(rpe_var), intent(in) :: mean, stdev
            type(rpe_var) :: u, v, randn
            real(dp) :: rand(2)

            call random_number(rand)

            ! Box-Muller method
            u = (-rpe_literal(2.0_dp) * log(rand(1))) ** rpe_literal(0.5_dp)
            v = rpe_literal(6.28318530718_dp) * rand(2)
            randn = mean + stdev * u * sin(v)
        end function randn

        !> @brief
        !> Seeds RNG from system clock.
        subroutine time_seed()
            integer :: i, n, clock
            integer, allocatable :: seed(:)

            call random_seed(size = n)
            allocate(seed(n))

            call system_clock(count=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call random_seed(put = seed)

            deallocate(seed)
        end subroutine time_seed
end module phy_sppt
