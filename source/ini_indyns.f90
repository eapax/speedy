subroutine indyns()
    ! subroutine indyns
    !
    ! Purpose : set time-stepping constants and initialize coefficients
    !           and spectral operators for model dynamics
    ! Initialized common blocks: mod_dyncon1
    !                            spectral (through routine parmtr)
    !

    use mod_tsteps, only: nsteps, alph
    use mod_dyncon0
    use mod_dyncon1
    use mod_atparam
    use mod_hdifcon, only: dmp, dmpd, dmps, tcorv, qcorv
    use spectral, only: sia, parmtr
    use mod_prec, only: dp

    implicit none

    integer :: j, k, jj, npowhd
    real(dp) :: elap, elapn, hdifd, hdiff, hdifs, qexp, rad1, rgam, rlap, twn

    ! 1. Definition of constants
    if (mod(nsteps,2)/=0) stop ' Invalid no. of time steps'

    ! alph = 0 ---- forward step for gravity wave terms
    ! alph = 1 ---- backward implicit -----------------
    ! alph = 0.5 -- centered implicit -----------------
    alph = 0.5_dp

    ! Power of Laplacian in horizontal diffusion
    npowhd = 4

    ! 2. Definition of model levels

    ! 2.1 Half (vertical velocity) levels
    if (kx==5) then
        hsg(:6) = (/ 0.000_dp, 0.150_dp, 0.350_dp, 0.650_dp, 0.900_dp, &
                1.000_dp /)
    else if (kx==7) then
        hsg(:8) = (/ 0.020_dp, 0.140_dp, 0.260_dp, 0.420_dp, 0.600_dp, &
                0.770_dp, 0.900_dp, 1.000_dp /)
    else if (kx==8) then
        hsg(:9) = (/ 0.000_dp, 0.050_dp, 0.140_dp, 0.260_dp, 0.420_dp, &
                0.600_dp, 0.770_dp, 0.900_dp, 1.000_dp /)
    end if

    do k = 1, kxp
        print *, ' Model half-level (*1000)', k, nint(HSG(k)*1000)
    end do

    ! 2.2 Layer thicknesses and full (u,v,T) levels
    do k = 1, kx
        dhs(k) = hsg(k+1)-hsg(k)
        fsg(k) = 0.5_dp*(hsg(k+1)+hsg(k))
    end do

    do k = 1, kx
        print *, ' Model full-level (*1000)', k, nint(FSG(k)*1000)
    end do

    ! 2.3 Additional functions of sigma
    do k = 1, kx
        dhsr(k) = 0.5_dp/dhs(k)
        fsgr(k) = akap/(2.0_dp*fsg(k))
    end do

    ! 3. Horizontal functions and spectral operators

    ! 3.1 Initialization of spectral operators
    call parmtr(1.0_dp)

    ! 3.2 Latitudes and functions of latitude
    !     NB: J=1 is Southernmost point!
    do j = 1, iy
        jj = il + 1 - j
        rad1 = asin(sia(j))
        radang(j)  = -rad1
        radang(jj) =  rad1
        coriol(j)  = -2.0_dp*omega*sia(j)*3600.0_dp
        coriol(jj) =  2.0_dp*omega*sia(j)*3600.0_dp
    end do

    ! 4. Coefficients to compute geopotential
    do k = 1, kx
      xgeop1(k) = rgas*log(hsg(k+1)/fsg(k))
      if (k/=kx) xgeop2(k+1) = rgas*log(fsg(k+1)/hsg(k+1))
    end do

    ! 5. Coefficients for horizontal diffusion

    ! 5.1 Spectral damping coefficients
    hdiff = 1.0_dp/(thd *3600.0_dp)
    hdifd = 1.0_dp/(thdd*3600.0_dp)
    hdifs = 1.0_dp/(thds*3600.0_dp)
    rlap  = 1.0_dp/float(mtrun*(mtrun+1))

    do j = 1, nx
        do k = 1, mx
            twn = float(isc*(k-1)+j-1)
            elap = (twn*(twn+1.0_dp)*rlap)
            elapn = elap**npowhd
            dmp(k,j)  = hdiff*elapn
            dmpd(k,j) = hdifd*elapn
            dmps(k,j) = hdifs*elap
        end do
    end do

    ! 5.2 Orographic correction terms for temperature and humidity
    !     (vertical component)
    rgam = rgas*gamma/(1000.0_dp*grav)
    qexp = hscale/hshum

    tcorv(1)=0.0_dp
    qcorv(1)=0.0_dp
    qcorv(2)=0.0_dp

    do k = 2, kx
        tcorv(k) = fsg(k)**rgam
        if (k>2) qcorv(k) = fsg(k)**qexp
        print *, ' temp/hum correction at level ', k, tcorv(k), qcorv(k)
    end do
end subroutine indyns
