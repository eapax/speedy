subroutine geop(jj)
    ! subroutine geop (jj)
    !
    ! Purpose : compute spectral geopotential from spectral temperature T
    !           and spectral topography PHIS, as in GFDL Climate Group GCM
    ! Input :   jj = time level index (1 or 2)
    ! Modified common blocks : DYNSP2

    use mod_atparam
    use mod_dynvar
    use mod_dyncon1, only: xgeop1, xgeop2, hsg, fsg
    use rp_emulator
    use mod_prec, only: set_precision

    implicit none

    integer, intent(in) :: jj
    integer :: k, m, n
    type(rpe_var) :: corf

    ! 1. Bottom layer (integration over half a layer)
    do n=1,nx
        do m=1,mx
            call set_precision(m, n)
            phi(m,n,kx) = phis(m,n) + xgeop1(kx) * t(m,n,kx,jj)
        end do
    end do

    ! 2. Other layers (integration two half-layers)
    do k = kx-1,1,-1
        do n=1,nx
            do m=1,mx
                call set_precision(m, n)
                phi(m,n,k) = phi(m,n,k+1) + xgeop2(k+1)*t(m,n,k+1,jj)&
                             & + xgeop1(k)*t(m,n,k,jj)
            end do
        end do
    end do

    ! 3. lapse-rate correction in the free troposphere
    do k = 2,kx-1
        do n=1,nx
            call set_precision(1, n)
            corf=xgeop1(k)*rpe_literal(0.5)*log(hsg(k+1)/fsg(k))/log(fsg(k+1)/fsg(k-1))
            phi(1,n,k) = phi(1,n,k) + corf*(t(1,n,k+1,jj) - t(1,n,k-1,jj))
        end do
    end do
    call set_precision(0, 0)
end
