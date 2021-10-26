subroutine forint(ngp,imon,fmon,for12,for1)
    ! Aux. routine FORINT : linear interpolation of monthly-mean forcing

    use rp_emulator
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: ngp, imon
    type(rpe_var), intent(in) :: fmon, for12(ngp,*)
    type(rpe_var), intent(inout) :: for1(ngp)
    integer :: imon2
    type(rpe_var) :: wmon, half

    half = 0.5_dp

    if (fmon<=half) then
        imon2 = imon-1
        if (imon==1) imon2 = 12
        wmon = half-fmon
    else
        imon2 = imon+1
        if (imon==12) imon2 = 1
        wmon = fmon-half
    end if

    for1 = for12(:,imon) + wmon*(for12(:,imon2) - for12(:,imon))
end subroutine forint

subroutine forin5(ngp,imon,fmon,for12,for1)
    ! Aux. routine FORIN5 : non-linear, mean-conserving interpolation
    !                       of monthly-mean forcing fields

    use rp_emulator
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: ngp, imon
    type(rpe_var), intent(in) :: fmon, for12(ngp,12)
    type(rpe_var), intent(inout) :: for1(ngp)
    integer :: im1, im2, ip1, ip2
    type(rpe_var) :: c0, t0, t1, t2, wm1, wm2, w0, wp1, wp2, one

    one = 1.0_dp

    im2 = imon-2
    im1 = imon-1
    ip1 = imon+1
    ip2 = imon+2

    if (im2<1)  im2 = im2+12
    if (im1<1)  im1 = im1+12
    if (ip1>12) ip1 = ip1-12
    if (ip2>12) ip2 = ip2-12

    c0 = one !/rpe_literal(12.0_dp)
    t0 = c0*fmon
    t1 = c0*(one-fmon)
    t2 = rpe_literal(0.25_dp)*fmon*(one-fmon)

    wm2 =        -t1   +t2
    wm1 =  -c0 +8*t1*rpe_literal(12.0_dp) -6*t2*rpe_literal(12.0_dp)
    w0  = 7*c0  +10*t2*rpe_literal(12.0_dp)
    wp1 =  -c0 +8*t0*rpe_literal(12.0_dp) -6*t2*rpe_literal(12.0_dp)
    wp2 =        -t0   +t2

    print *, 'forin5 weights = ', wm2,wm1,w0,wp1,wp2
    print *, 'wm2 deep dive:   ', wm2, -t1,t2, -t1+t2
    print *, '-t1 deep dive:   ', -t1,c0,one,fmon, c0*(one-fmon)
    print *, 'c0 deep dive:    ', c0, one,rpe_literal(12.0_dp),one/rpe_literal(12.0_dp),rpe_literal(1.0_dp)/rpe_literal(12.0_dp)  

    for1 = rpe_literal(12.0_dp)*wm2*for12(:,im2) + wm1*for12(:,im1) + w0*for12(:,imon) +&
        & wp1*for12(:,ip1) + rpe_literal(12.0_dp)*wp2*for12(:,ip2)

    for1 = for1/rpe_literal(12.0_dp)
end subroutine forin5
