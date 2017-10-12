!> @author
!> Sam Hatfield, AOPP, University of Oxford
!> @brief
!> A module for computing SPPT patterns to be used as multiplicative noise applied to physical tendencies
!> Stochastically Perturbed Parametrization Tendencies (SPPT) is a parametrization of model error.
!> See ECMWF Tech. Memo. #598 (Palmer et al. 2009)
module mod_sppt
    use mod_atparam
    use mod_tsteps, only: nsteps
    use mod_dyncon1, only: rearth
    use mod_spectral, only: el2

    implicit none

    private
    public mu, gen_sppt

    ! Array for tapering value of SPPT in the different layers of the atmosphere
    ! A value of 1 means the tendency is not tapered at that level
    real :: mu(kx) = (/ 1, 1, 1, 1, 1, 1, 1, 1 /)

    logical :: first = .true.

    ! Number of correlation scales for SPPT perturbations
    integer, parameter :: nscales = 3

    ! Decorrelation time of SPPT perturbation (in hours)
    real, dimension(nscales), parameter :: &
            time_decorr = (/ 3.0, 72.0, 720.0 /)

    ! Correlation length scale of SPPT perturbation (in metres)
    real, dimension(nscales), parameter :: &
            len_decorr = (/ 500000.0, 1000000.0, 2000000.0 /)

    ! Standard deviation of SPPT perturbation (in grid point space)
    real, dimension(nscales), parameter :: &
            stddev = (/ 0.52,  0.18, 0.06 /)

    ! Time autocorrelation of spectral AR(1) signals
    real :: phi(nscales)

    ! Total wavenumber-wise standard deviation of spectral signals
    real :: sigma(mx, nx, nscales)

    ! Perturbations in spectral space
    complex :: sppt_spec(mx, nx, nscales)

    contains
        !> @brief
        !> Generate grid point space SPPT pattern
        !> distribution.
        !> @return sppt_grid the generated grid point pattern
        function gen_sppt() result(sppt_grid_out)
            integer :: m, n, k
            complex :: sppt_spec_total(mx, nx)
            real :: sppt_grid(ix, il), sppt_grid_out(ix*il, kx)
            complex :: eta(mx, nx, nscales)
            real :: randreal, randimag

            ! Seed RNG if first use of SPPT
            if (first) then
                call time_seed()
                call init_sppt_parameters()
            end if

            ! Generate Gaussian noise
            do k=1, nscales
                do n=1, nx
                    do m=1, mx
                        randreal = randn(0.0, 1.0)
                        randimag = randn(0.0, 1.0)

                        ! Clip noise to +- 10 standard deviations
                        eta(m,n,k) = cmplx(&
                            & min(10.0, abs(randreal)) * sign(1.0,randreal),&
                            & min(10.0, abs(randimag)) * sign(1.0,randimag))
                    end do
                end do
            end do

            ! If first timestep
            if (first) then
                ! First AR(1) step
                do n=1, nscales
                    sppt_spec(:,:,n) = (1 - phi(n)**2)**(-0.5)
                end do
                sppt_spec = sppt_spec * sigma * eta
                first = .false.
            else
                ! Subsequent AR(1) steps
                do n=1, nscales
                    sppt_spec(:,:,n) = phi(n)*sppt_spec(:,:,n)
                end do
                sppt_spec = sppt_spec + sigma*eta
            end if

            ! Sum SPPT perturbations over correlation scales
            do n=1, nscales
                sppt_spec_total = sppt_spec_total + sppt_spec(:,:,n)
            end do

            ! Convert to grid point space
            call grid(sppt_spec_total, sppt_grid, 1)
            do k=1,kx
                ! SPPT perturbations uniform in height
                sppt_grid_out(:, k) = reshape(sppt_grid, (/ix*il/))
            end do

            ! Clip to +/- 1.0
            sppt_grid_out = min(1.0, abs(sppt_grid_out)) * sign(1.0,sppt_grid_out)
        end function

        !> @brief
        !> Calculates the phi and sigma parameters in the SPPT calculations
        !> from the defined correlation scales and standard deviations
        subroutine init_sppt_parameters()
            integer :: n
            real :: f0(nscales)

            ! Calculate time autocorrelation factor as a function of timestep
            phi = exp(-(24/real(nsteps))/time_decorr)

            ! Generate spatial amplitude pattern
            f0 = sum((/ ((2*n+1)*exp(-0.5*(len_decorr/rearth)**2*n*(n+1)),n=1,ntrun) /))
            f0 = sqrt((stddev**2*(1-phi**2))/f0)

            do n=1, nscales
                sigma(:,:,n) = f0(n) * exp(-0.25*len_decorr(n)**2 * el2)
            end do

        end subroutine

        !> @brief
        !> Generates a random number drawn for the specified normal
        !> distribution.
        !> @param mean the mean of the distribution to draw from
        !> @param stdev the standard deviation of the distribution to draw from
        !> @return randn the generated random number
        function randn(mean, stdev)
            real, intent(in) :: mean, stdev
            real :: u, v, randn
            real :: rand(2)

            call random_number(rand)

            ! Box-Muller method
            u = (-2.0 * log(rand(1))) ** 0.5
            v =   6.28318530718 * rand(2)
            randn = mean + stdev * u * sin(v)
        end function

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
        end subroutine
end module
