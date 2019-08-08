subroutine ini_atm()
    !   subroutine ini_atm (cexp)
    !
    !   purpose : call initialization routines for all model common blocks

    use mod_atparam
    use mod_dyncon1, only: grav, hsg, fsg, radang
    use mod_fluxes, only: ini_fluxes
    use mod_prec, only: dp

    implicit none

    real(dp) :: ppl(kx)            ! post-processing levels (hpa/1000)
    integer :: k

    ! 1. initialize ffts
    call inifft()

    ! 2. initialize dynamical constants and operators
    call indyns()

    ! 3. set post-processing levels
    do k = 1, kx
        ppl(k) = prlev(fsg(k))
    end do

    ! 4. initialize forcing fields (boundary cond.)
    call inbcon(grav,radang)

    ! 5. initialize model variables
    call invars()

    ! 6. initialize constants for physical parametrization
    call inphys(hsg, ppl, radang)

    ! 7. initialize time-mean arrays for surface fluxes and output fields
    call ini_fluxes()

    contains
        function prlev(siglev)
            ! function prlev (siglev)
            ! purpose : select the closest standard pressure level for post-proc.
            ! input :   siglev = sigma level
            use mod_prec, only: dp

            implicit none

            real(dp), intent(in) :: siglev
            real(dp) :: plev(14)
            real(dp) :: prlev, dif, adif
            integer :: k

            plev = (/ 0.925_dp, 0.850_dp, 0.775_dp, 0.700_dp, 0.600_dp, &
                    0.500_dp, 0.400_dp, 0.300_dp, 0.250_dp, 0.200_dp, &
                    0.150_dp, 0.100_dp, 0.050_dp, 0.030_dp /)

            dif = 1.0_dp - siglev

            prlev = 1.0_dp

            do k = 1, 14
                adif = abs(plev(k) - siglev)
                if (adif<=dif) then
                    dif = adif
                    prlev = plev(k)
                end if
            end do
        end function prlev
end subroutine ini_atm
