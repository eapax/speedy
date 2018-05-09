module spectral

    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    private
    public setup_spectral, gaussl, parmtr, lap, invlap, grad, uvspec, grid, &
            spec, vdspec
    public el2, sia, cosg, cosgr
    public cpol, nsh2, wt, trfilt

    ! Initial. in parmtr
    type(rpe_var), dimension(:,:), allocatable :: el2, elm2, el4, trfilt
    integer, allocatable :: l2(:,:), ll(:,:), mm(:), nsh2(:)

    ! Initial. in parmtr
    type(rpe_var), dimension(:), allocatable :: sia, coa, wt, wght
    type(rpe_var), dimension(:), allocatable :: cosg, cosgr, cosgr2

    ! Initial. in parmtr
    type(rpe_var), allocatable :: gradx(:), gradym(:,:), gradyp(:,:)

    ! Initial. in parmtr
    type(rpe_var) :: sqrhlf
    type(rpe_var), allocatable :: consq(:), epsi(:,:), repsi(:,:), emm(:), ell(:,:)

    ! Initial. in parmtr
    type(rpe_var), allocatable :: poly(:,:)

    ! Initial. in parmtr
    type(rpe_var), allocatable :: cpol(:,:,:)

    ! Initial. in parmtr
    type(rpe_var), dimension(:,:), allocatable :: uvdx, uvdym, uvdyp

    ! Initial. in parmtr
    type(rpe_var), dimension(:,:), allocatable :: vddym, vddyp

    contains
        subroutine setup_spectral()
            allocate(el2(mx, nx))
            allocate(elm2(mx, nx))
            allocate(el4(mx, nx))
            allocate(trfilt(mx, nx))
            allocate(l2(mx,nx))
            allocate(ll(mx,nx))
            allocate(mm(mx))
            allocate(nsh2(nx))
            allocate(sia(iy))
            allocate(coa(iy))
            allocate(wt(iy))
            allocate(wght(iy))
            allocate(cosg(il))
            allocate(cosgr(il))
            allocate(cosgr2(il))
            allocate(gradx(mx))
            allocate(gradym(mx, nx))
            allocate(gradyp(mx, nx))
            allocate(consq(mxp))
            allocate(epsi(mxp,nxp))
            allocate(repsi(mxp,nxp))
            allocate(emm(mxp))
            allocate(ell(mxp,nxp))
            allocate(poly(mx,nx))
            allocate(cpol(mx2,nx,iy))
            allocate(uvdx(mx, nx))
            allocate(uvdym(mx, nx))
            allocate(uvdyp(mx, nx))
            allocate(vddym(mx, nx))
            allocate(vddyp(mx, nx))
        end subroutine setup_spectral

!******************************************************************
subroutine gaussl(x,w,m)
    !   a slightly modified version of a program in Numerical Recipes 
    !       (Cambridge Univ. Press, 1989)
    !   input:
    !      m    = number of gaussian latitudes between pole and equator
    !   output:
    !      x(m) = sin(gaussian latitude) 
    !      w(m) = weights in gaussian quadrature (sum should equal 1.0)

    type(rpe_var), intent(inout) :: x(m),w(m)
    integer, intent(in) :: m
    real(dp) :: z,z1,p1,p2,p3,pp
    real(dp), parameter :: eps=3.d-14
    integer :: n, j, i

    n = 2*m

    z1 = 2.0_dp

    do i=1,m
        z=cos(3.141592654_dp*(i-.25_dp)/(n+.5_dp))
        do while (abs(z-z1).gt.eps)
            p1=1._dp
            p2=0._dp
    
            do j=1,n
              p3=p2
              p2=p1
              p1=((2._dp*j-1._dp)*z*p2-(j-1._dp)*p3)/j
            end do
    
            pp=n*(z*p1-p2)/(z*z-1._dp)
            z1=z
            z=z1-p1/pp
        end do

        x(i)=z
        w(i)=2.0_dp/((1.0_dp-z*z)*pp*pp)
    end do
end
!****************************************************************
subroutine parmtr(a)

    type(rpe_var), intent(in) :: a
    type(rpe_var) :: am1, am2, cosqr, el1, ell2, emm2

    integer :: j, jj, m, m1, m2, n

    ! initializes Legendre transforms and constants used for other
    ! subroutines that manipulate spherical harmonics
    !
    ! input:  A = radius of the sphere
    ! first compute Gaussian latitudes and weights at the IY points from 
    !     pole to equator
    ! SIA(IY) is sin of latitude, WT(IY) are Gaussian weights for quadratures,
    !   saved in mod_spectral
    call gaussl(sia,wt,iy)
    am1 = 1.0_dp/a
    am2=  1.0_dp/(a*a)

    ! COA(IY) = cos(lat); WGHT needed for transforms, 
    !           saved in mod_spectral
    do j=1,iy
        cosqr = 1.0_dp-sia(j)**2
        coa(j)=sqrt(cosqr)
        wght(j)=wt(j)/(a*cosqr)
    end do

    ! expand cosine and its reciprocal to cover both hemispheres, 
    !    saved in mod_spectral
    do j=1,iy
        jj=il+1-j
        cosg(j)=coa(j)
        cosg(jj)=coa(j)
        cosgr(j)=1.0_dp/coa(j)
        cosgr(jj)=1.0_dp/coa(j)
        cosgr2(j)=1.0_dp/(coa(j)*coa(j))
        cosgr2(jj)=1.0_dp/(coa(j)*coa(j))
    end do

    !  MM = zonal wavenumber = m
    !     ISC=3 implies that only wavenumber 0,3,6,9,etc are included in model
    !  LL = total wavenumber of spherical harmonic = l
    !  L2 = l*(l+1)
    !  EL2 = l*(l+1)/(a**2)
    !  EL4 = EL2*EL2 ; for biharmonic diffusion
    !  ELM2 = 1./EL2
    !  TRFILT used to filter out "non-triangular" part of rhomboidal truncation
    !   saved in mod_spectral
    do n=1,nx
        nsh2(n)=0
        do m=1,mx
            mm(m)=isc*(m-1)
            ll(m,n)=mm(m)+n-1
            l2(m,n)=ll(m,n)*(ll(m,n)+1)
            el2(m,n)=float(l2(m,n))*am2
            el4(m,n)=el2(m,n)*el2(m,n)
            if (ll(m,n).le.ntrun1.or.ix.ne.4*iy) nsh2(n)=nsh2(n)+2
            if (ll(m,n).le.ntrun) then
              trfilt(m,n)=1.0_dp
            else
              trfilt(m,n)=0.0_dp
            end if
        end do
    end do

    elm2(1,1)=0.0_dp
    do m=2,mx
        do n=1,nx
            elm2(m,n)=1.0_dp/el2(m,n)
        end do
    end do

    do n=2,nx
        elm2(1,n)=1.0_dp/el2(1,n)
    end do

    ! quantities needed to generate and differentiate Legendre polynomials
    ! all m values up to MXP = ISC*MTRUN+1 are needed by recursion relation 
    ! saved in mod_spectral
    do m=1,mxp
        do n=1,nxp
            emm(m)=float(m-1)
            ell(m,n)=float(n+m-2)
            emm2=emm(m)**2
            ell2=ell(m,n)**2
            if(n.eq.nxp) then
              epsi(m,n)=0.0_dp
            else if(n.eq.1.and.m.eq.1) then
              epsi(m,n)=0.0_dp
            else
              epsi(m,n)=sqrt((ell2-emm2)/(4.0_dp*ell2-1.0_dp))
            end if
            repsi(m,n)=0.0_dp
            if(epsi(m,n).gt.0.0_dp) repsi(m,n)=1.0_dp/epsi(m,n)
        end do
    end do

    sqrhlf=sqrt(0.5_dp)
    do m=2,mxp
        consq(m) = sqrt(0.5_dp*(2.0_dp*emm(m)+1.0_dp)/emm(m))
    end do

    ! quantities required by subroutines GRAD, UVSPEC, and VDS
    ! saved in mod_spectral
    do m=1,mx
        do n=1,nx
            m1=mm(m)
            m2=m1+1
            el1=float(ll(m,n))
            if(n.eq.1) then
                gradx(m)=float(m1)/a
                uvdx(m,1)=-a/float(m1+1)
                uvdym(m,1)=0.0_dp
                vddym(m,1)=0.0_dp
            else
                uvdx(m,n)=-a*float(m1)/(el1*(el1+1))
                gradym(m,n)=(el1-1.0_dp)*epsi(m2,n)/a
                uvdym(m,n)=-a*epsi(m2,n)/el1
                vddym(m,n)=(el1+1)*epsi(m2,n)/a
            end if
            gradyp(m,n)=(el1+2.0_dp)*epsi(m2,n+1)/a
            uvdyp(m,n)=-a*epsi(m2,n+1)/(el1+1.0_dp)
            vddyp(m,n)=el1*epsi(m2,n+1)/a
        end do
    end do

    !  generate associated Legendre polynomial
    !  LGNDRE computes the polynomials at a particular latitiude, POLY(MX,NX), and stores
    !  them in mod_spectral
    !  polynomials and 'clones' stored in mod_spectral
    do j=1,iy
        call lgndre(j)
        do n=1,nx
            do m=1,mx
                m1=2*m-1
                m2=2*m
                cpol(m1,n,j)=poly(m,n)
                cpol(m2,n,j)=poly(m,n)
            end do
        end do
    end do
end
!****************************************************************
subroutine lgndre(j)
    ! follows Leith Holloways code

    integer, intent(in) :: j
    real(dp), parameter :: small = 1.0d-30

    integer :: m, n, mm2
    type(rpe_var) :: alp(mxp,nx), x, y
    y = coa(j)
    x = sia(j)

    ! start recursion with N=1 (M=L) diagonal 
    alp(1,1) = sqrhlf
    do m=2,mxp
        alp(m,1) = consq(m)*y*alp(m-1,1)
    end do
  
    ! continue with other elements
    do m=1,mxp
        alp(m,2)=(x*alp(m,1))*repsi(m,2)
    end do

    do n=3,nx
        do m=1,mxp
          alp(m,n)=(x*alp(m,n-1)-epsi(m,n-1)*alp(m,n-2))*repsi(m,n)
        end do
    end do

    ! zero polynomials with absolute values smaller than 10**(-30)
    do n=1,nx
        do m=1,mxp
            if(abs(alp(m,n)) .le. small) alp(m,n)=0.0_dp
        end do
    end do

    ! pick off the required polynomials
    do n=1,nx
        do m=1,mx
            mm2=isc*(m-1)+1
            poly(m,n)=alp(mm2,n)
        end do
    end do
end
!***************************************************************
subroutine lap(strm,vorm)
    type(rpe_complex_var), intent(in) :: strm(mx,nx)
    type(rpe_complex_var), intent(inout) :: vorm(mx,nx)

    vorm = -strm * el2
end
!*******************************************************************
subroutine invlap(vorm,strm)
    type(rpe_complex_var), intent(in) :: vorm(mx,nx)
    type(rpe_complex_var), intent(inout) :: strm(mx,nx)

    strm = -vorm * elm2
end
!*********************************************************************
subroutine grad(psi,psdx,psdy)
    type(rpe_complex_var), dimension(mx,nx), intent(inout) :: psi
    type(rpe_complex_var), dimension(mx,nx), intent(inout) :: psdx, psdy

    integer :: n

    do n=1,nx
        psdx(:,n) = CMPLX(-gradx*REAL(AIMAG(psi(:,n))), &
                           gradx*REAL(REAL (psi(:,n))))
    end do

    psdy(:,1) = gradyp(:,1)*psi(:,2)
    psdy(:,nx) = -gradym(:,nx)*psi(:,ntrun1)


    do n=2,ntrun1
        psdy(:,n) = -gradym(:,n)*psi(:,n-1) + gradyp(:,n)*psi(:,n+1)
    end do
end
!******************************************************************
subroutine vds(ucosm,vcosm,vorm,divm)
    type(rpe_complex_var), dimension(mx,nx) :: ucosm, vcosm
    type(rpe_complex_var), dimension(mx,nx), intent(inout) :: vorm, divm
    type(rpe_complex_var), dimension(mx,nx) :: zc, zp
    
    integer :: n

    do n=1,nx
        zp(:,n) = CMPLX(-gradx*REAL(AIMAG(ucosm(:,n))), &
                        gradx*REAL(REAL (ucosm(:,n))))
        zc(:,n) = CMPLX(-gradx*REAL(AIMAG(vcosm(:,n))), &
                        gradx*REAL(REAL (vcosm(:,n))))
    end do

    vorm(:,1) = zc(:,1) - vddyp(:,1)*ucosm(:,2)
    vorm(:,nx) = vddym(:,nx)*ucosm(:,ntrun1)
    divm(:,1) = zp(:,1) + vddyp(:,1)*vcosm(:,2)
    divm(:,nx) = -vddym(:,nx)*vcosm(:,ntrun1)

    do n=2,ntrun1
        vorm(:,n)= vddym(:,n)*ucosm(:,n-1) - vddyp(:,n)*ucosm(:,n+1) + zc(:,n)
        divm(:,n)=-vddym(:,n)*vcosm(:,n-1) + vddyp(:,n)*vcosm(:,n+1) + zp(:,n)
    end do
end
!******************************************************************
subroutine uvspec(vorm, divm, um, vm)
    ! Calculate u and v in grid-point space from vorticity and divergence in
    ! spectral space
    type(rpe_complex_var), dimension(mx,nx), intent(in) :: vorm, divm
    type(rpe_var), dimension(ix,il), intent(out) :: um, vm
    type(rpe_complex_var), dimension(mx,nx) :: ucosm, vcosm
    type(rpe_complex_var), dimension(mx,nx) :: zc, zp

    integer :: n

    zp = CMPLX(-uvdx*REAL(AIMAG(vorm)), uvdx*REAL(REAL(vorm)))
    zc = CMPLX(-uvdx*REAL(AIMAG(divm)), uvdx*REAL(REAL(divm)))

    ucosm(:,1) = zc(:,1) - uvdyp(:,1)*vorm(:,2)
    ucosm(:,nx) = uvdym(:,nx)*vorm(:,ntrun1)
    vcosm(:,1) = zp(:,1) + uvdyp(:,1)*divm(:,2)
    vcosm(:,nx) = -uvdym(:,nx)*divm(:,ntrun1)

    do n=2,ntrun1
        vcosm(:,n) =-uvdym(:,n)*divm(:,n-1) + uvdyp(:,n)*divm(:,n+1) + zp(:,n)
        ucosm(:,n) = uvdym(:,n)*vorm(:,n-1) - uvdyp(:,n)*vorm(:,n+1) + zc(:,n)
    end do

    call grid(ucosm, um, 2)
    call grid(vcosm, vm, 2)
end
!*******************************************************************
subroutine grid(vorm,vorg,kcos)
    type(rpe_var), intent(out) :: vorg(ix,il)
    type(rpe_complex_var), intent(in) :: vorm(mx,nx)
    integer, intent(in) :: kcos
    type(rpe_complex_var) :: varm(mx,il)

    call gridy(vorm,varm)
    call gridx(varm,vorg,kcos)
end
!*********************************************************************
subroutine spec(vorg,vorm)
    type(rpe_var), intent(in) :: vorg(ix,il)
    type(rpe_complex_var), intent(out) :: vorm(mx,nx)
    type(rpe_complex_var) :: varm(mx,il)

    call specx(vorg,varm)
    call specy(varm,vorm)
end
!*********************************************************************
subroutine vdspec(ug,vg,vorm,divm,kcos)

    type(rpe_complex_var), intent(in) :: ug(ix,il), vg(ix,il)
    type(rpe_complex_var), intent(out) :: vorm(mx,nx), divm(mx,nx)
    integer, intent(in) :: kcos
    integer :: i, j
    type(rpe_var) :: ug1(ix,il), vg1(ix,il)
    type(rpe_complex_var) :: um(mx,il), vm(mx,il), dumc1(mx,nx), dumc2(mx,nx)

    if (kcos.eq.2) then
        do j=1,il
            do i=1,ix
                ug1(i,j)=ug(i,j)*cosgr(j)
                vg1(i,j)=vg(i,j)*cosgr(j)
            end do
        end do
    else
        do j=1,il
            do i=1,ix
                ug1(i,j)=ug(i,j)*cosgr2(j)
                vg1(i,j)=vg(i,j)*cosgr2(j)
            end do
        end do
    end if

    call specx(ug1,um)  
    call specx(vg1,vm)
    call specy(um,dumc1)
    call specy(vm,dumc2)
    call vds(dumc1,dumc2,vorm,divm)
end

end module spectral
!*********************************************************************
subroutine gridy(v,varm)
    use mod_atparam
    use spectral, only: cpol, nsh2

    implicit none

    type(rpe_var), intent(in) :: v(mx2,nx)
    type(rpe_var), intent(inout) :: varm(mx2,il)
    type(rpe_var) :: vm1(mx2),vm2(mx2)

    integer :: j, j1, m, n

    do j=1,iy
        j1=il+1-j

        do m=1,mx2
            vm1(m)=0.0_dp
            vm2(m)=0.0_dp
        end do

        do n=1,nx,2
            !do m=1,mx2
            do m=1,nsh2(n)
                vm1(m)=vm1(m)+v(m,n)*cpol(m,n,j)
            end do
        end do

        do n=2,nx,2
            !do m=1,mx2
            do m=1,nsh2(n)
                vm2(m)=vm2(m)+v(m,n)*cpol(m,n,j)
            end do
        end do

        do m=1,mx2
            varm(m,j1)=vm1(m)+vm2(m)
            varm(m,j) =vm1(m)-vm2(m)
        end do
    end do
end
!******************************************************************
subroutine specy(varm,vorm)
    use mod_atparam
    use spectral, only: wt, cpol, nsh2
    use rp_emulator
    use mod_prec

    implicit none

    type(rpe_var), intent(in) :: varm(mx2,il)
    type(rpe_var), intent(inout) :: vorm(mx2,nx)
    type(rpe_var) :: svarm(mx2,iy), dvarm(mx2,iy)

    integer :: j, j1, m, n

    vorm = 0.0_dp

    do j=1,iy
        j1=il+1-j
        do m=1,mx2
            svarm(m,j)=(varm(m,j1)+varm(m,j))*wt(j)
            dvarm(m,j)=(varm(m,j1)-varm(m,j))*wt(j)
        end do
    end do

    do j=1,iy
        j1=il+1-j

        do n=1,ntrun1,2
            !do m=1,mx2
            do m=1,nsh2(n)
                vorm(m,n) = vorm(m,n)+cpol(m,n,j)*svarm(m,j)
            end do
        end do

        do n=2,ntrun1,2
            !do m=1,mx2
            do m=1,nsh2(n)
                vorm(m,n) = vorm(m,n)+cpol(m,n,j)*dvarm(m,j)
            end do
        end do
    end do
end
!******************************************************************
subroutine trunct(vor)
    use mod_atparam
    use spectral, only: trfilt
    use rp_emulator

    implicit none

    type(rpe_complex_var), intent(inout) :: vor(mx,nx)

    vor = vor * trfilt
end
