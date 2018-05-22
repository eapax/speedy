subroutine ini_atm()
    !   subroutine ini_atm (cexp)
    !
    !   purpose : call initialization routines for all model common blocks 

    use mod_atparam
    use mod_dyncon1, only: grav, hsg, fsg, radang
    use phy_sppt, only: ini_sppt

    implicit none

    real :: ppl(kx)            ! post-processing levels (hpa/1000)
    integer :: iitest = 1, k

    ! 1. initialize ffts
    if (iitest == 1) print *, 'calling inifft'
    call inifft()

    ! 2. initialize dynamical constants and operators
    if (iitest == 1) print *, 'calling indyns'
    call indyns()

    ! 3. set post-processing levels
    do k = 1, kx
        ppl(k) = prlev(fsg(k))
    end do

    ! 4. initialize constants for physical parametrization
    if (iitest == 1) print *, 'calling inphys'
    call inphys(hsg, ppl, radang)
    call ini_sppt()

    ! 5. initialize forcing fields (boundary cond.)
    if (iitest == 1) print *, 'calling inbcon'
    call inbcon(grav,radang)

    ! 6. initialize model variables
    if (iitest == 1) print *, 'calling invars'
    call invars()

    ! 7. initialize time-mean arrays for surface fluxes and output fields
    if (iitest == 1) print *, 'calling dmflux'
    call dmflux(0)

    contains
        function prlev(siglev)
            ! function prlev (siglev)
            ! purpose : select the closest standard pressure level for post-proc.
            ! input :   siglev = sigma level
            implicit none
        
            real, intent(in) :: siglev
            real :: plev(14) = (/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,&
                & 0.300, 0.250, 0.200, 0.150, 0.100, 0.050, 0.030 /)
            real :: prlev, dif, adif
            integer :: k
        
            dif = 1.0 - siglev
        
            prlev = 1.0
        
            do k = 1, 14
                adif = abs(plev(k) - siglev)
                if (adif <= dif) then
                    dif = adif
                    prlev = plev(k)
                end if
            end do
        end
end
