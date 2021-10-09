subroutine geop(jj)
    ! subroutine geop (jj)
    !
    ! Purpose : compute spectral geopotential from spectral temperature T
    !           and spectral topography PHIS, as in GFDL Climate Group GCM
    ! Input :   jj = time level index (1 or 2)
    ! Modified common blocks : DYNSP2

    use mod_atparam
    use mod_dynvar, only: t, phi, phis, tcopy
    use mod_dyncon1, only: xgeop1, xgeop2, hsg, fsg
    use humidity, only: zero_C
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: jj
    integer :: k
    type(rpe_var) :: corf

    !Create a copy of the initial state of the spectral temperature
!    tcopy = t

    ! Convert temperature to Kelvin
    t(1,1,:,:) = t(1,1,:,:) + cmplx(sqrt(2.0_dp)*zero_c, kind=dp)

    ! 1. Bottom layer (integration over half a layer)
    phi(:,:,kx) = phis + xgeop1(kx) * t(:,:,kx,jj)

    ! 2. Other layers (integration two half-layers)
    do k = kx-1,1,-1
        phi(:,:,k) = phi(:,:,k+1) + xgeop2(k+1)*t(:,:,k+1,jj) + &
                xgeop1(k)*t(:,:,k,jj)
    end do

    ! 3. lapse-rate correction in the free troposphere
    do k = 2,kx-1
        corf=xgeop1(k)*rpe_literal(0.5_dp)* &
                log(hsg(k+1)/fsg(k)) / log(fsg(k+1)/fsg(k-1))
        phi(1,:,k) = phi(1,:,k) + corf*(t(1,:,k+1,jj) - t(1,:,k-1,jj))
    end do

    ! Convert temperature to Celsius
    t(1,1,:,:) = t(1,1,:,:) - cmplx(sqrt(2.0_dp)*zero_c, kind=dp)
    
    !Rather than doing the previous line, just return the copy we created at the start
    !t = tcopy


end subroutine geop
