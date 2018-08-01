subroutine inphys(hsg,ppl,rlat)
    !
    ! subroutine inphys (hsg,ppl,rlat)
    !
    ! Purpose: Initialize mod_physvar
    ! Input :  hsg  : sigma at half levels
    !          ppl  : pressure levels for post-processing
    !          rlat : gaussian-grid latitudes

    use mod_atparam
    use mod_physcon
    use phy_convmf, only: ini_convmf
    use phy_lscond, only: ini_lscond
    use phy_vdifsc, only: ini_vdifsc
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    type(rpe_var), intent(in) :: hsg(0:kx)
    real(dp), intent(in) :: ppl(kx)
    type(rpe_var), intent(in) :: rlat(il)
    integer :: j, k

    ! 1.2 Functions of sigma and latitude
    sigh(0) = hsg(0)

    do k = 1, kx
        sig(k)  = 0.5_dp*(hsg(k)+hsg(k-1))
        sigl(k) = log(sig(k))
        sigh(k) = hsg(k)
        dsig(k) = hsg(k)-hsg(k-1)
        pout(k) = ppl(k)
        grdsig(k) = gg/(dsig(k)*p0)
        grdscp(k) = grdsig(k)/cp
    end do

    ! Weights for vertical interpolation at half-levels(1,kx) and surface
    ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
    ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
    ! Fsurf = Ffull(kx)+WVI(kx,2)*(Ffull(kx)-Ffull(kx-1))
    do k = 1, kx-1
        wvi(k,1) = 1.0_dp/(sigl(k+1)-sigl(k))
        wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
    end do

    wvi(kx,1) = 0.0_dp
    wvi(kx,2) = (log(0.99_dp)-sigl(kx))*wvi(kx-1,1)

    do j = 1, il
        slat(j) = sin(rlat(j))
        clat(j) = cos(rlat(j))
    end do

    ! Call setup routines for individual physics schemes
    call ini_convmf()
    call ini_lscond()
    call ini_vdifsc()
end subroutine inphys
