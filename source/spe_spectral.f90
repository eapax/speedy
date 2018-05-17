module spectral

    use mod_atparam

    implicit none

    private
    ! Public subroutines
    public setup_spectral, gaussl, parmtr, lap, invlap, grad, uvspec, grid, &
            spec, vdspec
    ! Public parameters
    public el2, sia, cosgr

    ! The variables below are declared public as they are used in the spectral
    ! subroutines at the bottom of this document that are outside the module
    ! due to issues with passing different types between subroutines but could
    ! potentially be fixed to private
    public cpol, nsh2, wt, trfilt

    ! Initial. in parmtr
    real, dimension(:,:), allocatable :: el2, elm2, trfilt
    integer, allocatable :: nsh2(:)

    ! Initial. in parmtr
    real, dimension(:), allocatable :: sia, wt
    real, dimension(:), allocatable :: cosgr, cosgr2

    ! Initial. in parmtr
    real, allocatable :: gradx(:), gradym(:,:), gradyp(:,:)

    ! Initial. in parmtr
    real, allocatable :: cpol(:,:,:)

    ! Initial. in parmtr
    real, dimension(:,:), allocatable :: uvdx, uvdym, uvdyp

    ! Initial. in parmtr
    real, dimension(:,:), allocatable :: vddym, vddyp

    contains
        subroutine setup_spectral()
            allocate(el2(mx, nx))
            allocate(elm2(mx, nx))
            allocate(trfilt(mx, nx))
            allocate(nsh2(nx))
            allocate(sia(iy))
            allocate(wt(iy))
            allocate(cosgr(il))
            allocate(cosgr2(il))
            allocate(gradx(mx))
            allocate(gradym(mx, nx))
            allocate(gradyp(mx, nx))
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

            real, intent(inout) :: x(m),w(m)
            integer, intent(in) :: m
            double precision :: z,z1,p1,p2,p3,pp
            double precision, parameter :: eps=3.d-14
            integer :: n, j, i

            n = 2*m

            z1 = 2.0

            do i=1,m
                z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
                do while (abs(z-z1).gt.eps)
                    p1=1.d0
                    p2=0.d0

                    do j=1,n
                      p3=p2
                      p2=p1
                      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
                    end do

                    pp=n*(z*p1-p2)/(z*z-1.d0)
                    z1=z
                    z=z1-p1/pp
                end do

                x(i)=z
                w(i)=2.d0/((1.d0-z*z)*pp*pp)
            end do
        end subroutine gaussl
        !****************************************************************
        subroutine parmtr(a)
            ! initializes Legendre transforms and constants used for other
            ! subroutines that manipulate spherical harmonics
            !
            ! input:  A = radius of the sphere
            ! first compute Gaussian latitudes and weights at the IY points from
            !     pole to equator
            ! SIA(IY) is sin of latitude, WT(IY) are Gaussian weights for
            !     quadratures
            real, intent(in) :: a
            real :: poly(mx, nx), coa(iy)
            integer :: l2(mx,nx), ll(mx,nx), mm(mx)
            real :: am1, am2, cosqr, el1, ell2, emm2, sqrhlf
            real :: consq(mxp), epsi(mxp,nxp), repsi(mxp,nxp), emm(mxp), &
                    ell(mxp,nxp)

            integer :: j, jj, m, m1, m2, n

            call gaussl(sia,wt,iy)
            am1 = 1./a
            am2=  1./(a*a)

            ! COA(IY) = cos(lat)
            do j=1,iy
                cosqr = 1.0-sia(j)**2
                coa(j)=sqrt(cosqr)
            end do

            ! expand cosine and its reciprocal to cover both hemispheres,
            !    saved in mod_spectral
            do j=1,iy
                jj=il+1-j
                cosgr(j)=1./coa(j)
                cosgr(jj)=1./coa(j)
                cosgr2(j)=1./(coa(j)*coa(j))
                cosgr2(jj)=1./(coa(j)*coa(j))
            end do

            !  MM = zonal wavenumber = m
            !    ISC=3 implies that only wavenumber 0,3,6,9,etc
            !    are included in model
            !  LL = total wavenumber of spherical harmonic = l
            !  L2 = l*(l+1)
            !  EL2 = l*(l+1)/(a**2)
            !  EL4 = EL2*EL2 ; for biharmonic diffusion
            !  ELM2 = 1./EL2
            !  TRFILT used to filter out "non-triangular" part of rhomboidal
            !  truncation
            do n=1,nx
                nsh2(n)=0
                do m=1,mx
                    mm(m)=isc*(m-1)
                    ll(m,n)=mm(m)+n-1
                    l2(m,n)=ll(m,n)*(ll(m,n)+1)
                    el2(m,n)=float(l2(m,n))*am2
                    if (ll(m,n).le.ntrun1.or.ix.ne.4*iy) nsh2(n)=nsh2(n)+2
                    if (ll(m,n).le.ntrun) then
                      trfilt(m,n)=1.
                    else
                      trfilt(m,n)=0.
                    end if
                end do
            end do

            elm2(1,1)=0.
            do m=2,mx
                do n=1,nx
                    elm2(m,n)=1./el2(m,n)
                end do
            end do

            do n=2,nx
                elm2(1,n)=1./el2(1,n)
            end do

            ! quantities needed to generate and differentiate Legendre
            ! polynomials. All m values up to MXP = ISC*MTRUN+1 are needed by
            ! recursion relation
            do m=1,mxp
                do n=1,nxp
                    emm(m)=float(m-1)
                    ell(m,n)=float(n+m-2)
                    emm2=emm(m)**2
                    ell2=ell(m,n)**2
                    if(n.eq.nxp) then
                      epsi(m,n)=0.0
                    else if(n.eq.1.and.m.eq.1) then
                      epsi(m,n)=0.0
                    else
                      epsi(m,n)=sqrt((ell2-emm2)/(4.*ell2-1.))
                    end if
                    repsi(m,n)=0.0
                    if(epsi(m,n).gt.0.) repsi(m,n)=1./epsi(m,n)
                end do
            end do

            do m=2,mxp
                consq(m) = sqrt(.5*(2.*emm(m)+1.)/emm(m))
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
                        uvdym(m,1)=0.0
                        vddym(m,1)=0.0
                    else
                        uvdx(m,n)=-a*float(m1)/(el1*(el1+1))
                        gradym(m,n)=(el1-1.)*epsi(m2,n)/a
                        uvdym(m,n)=-a*epsi(m2,n)/el1
                        vddym(m,n)=(el1+1)*epsi(m2,n)/a
                    end if
                    gradyp(m,n)=(el1+2.)*epsi(m2,n+1)/a
                    uvdyp(m,n)=-a*epsi(m2,n+1)/(el1+1.)
                    vddyp(m,n)=el1*epsi(m2,n+1)/a
                end do
            end do

            !  generate associated Legendre polynomial
            !  LGNDRE computes the polynomials at a particular latitiude
            !  POLY(MX,NX). These are then stored in cpol.
            sqrhlf=sqrt(.5)
            do j=1,iy
                call lgndre(j, sia(j), coa(j), poly, &
                        sqrhlf, consq, epsi, repsi)
                do n=1,nx
                    do m=1,mx
                        m1=2*m-1
                        m2=2*m
                        cpol(m1,n,j)=poly(m,n)
                        cpol(m2,n,j)=poly(m,n)
                    end do
                end do
            end do
        end subroutine parmtr
        !****************************************************************
        subroutine lgndre(j, x, y, poly, sqrhlf, consq, epsi, repsi)
            ! follows Leith Holloways code
            ! Sets the values of Legendre polynomial (poly)

            integer, intent(in) :: j
            real, intent(in) :: x, y
            real, intent(in) :: sqrhlf
            real, intent(in) :: consq(mxp), epsi(mxp,nxp), repsi(mxp,nxp)
            real, intent(out) :: poly(mx, nx)
            real, parameter :: small = 1.e-30

            integer :: m, n, mm2
            real :: alp(mxp,nx)

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
                    if(abs(alp(m,n)) .le. small) alp(m,n)=0.0
                end do
            end do

            ! pick off the required polynomials
            do n=1,nx
                do m=1,mx
                    mm2=isc*(m-1)+1
                    poly(m,n)=alp(mm2,n)
                end do
            end do
        end subroutine lgndre
        !***************************************************************
        subroutine lap(strm,vorm)
            ! Laplacian in spectral space

            complex, intent(in) :: strm(mx,nx)
            complex, intent(inout) :: vorm(mx,nx)

            vorm = -strm * el2
        end subroutine lap
        !*******************************************************************
        subroutine invlap(vorm,strm)
            ! Inverse Laplacian in spectral space

            complex, intent(in) :: vorm(mx,nx)
            complex, intent(inout) :: strm(mx,nx)

            strm = -vorm * elm2
        end subroutine invlap
        !*********************************************************************
        subroutine grad(psi,psdx,psdy)
            ! Gradient in spectral space

            complex, dimension(mx,nx), intent(inout) :: psi
            complex, dimension(mx,nx), intent(inout) :: psdx, psdy

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
        end subroutine grad
        !******************************************************************
        subroutine vds(ucosm,vcosm,vorm,divm)
            complex, dimension(mx,nx) :: ucosm, vcosm
            complex, dimension(mx,nx), intent(inout) :: vorm, divm
            complex, dimension(mx,nx) :: zc, zp

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
        end subroutine vds
        !******************************************************************
        subroutine uvspec(vorm, divm, um, vm)
            ! Calculate u and v in grid-point space from vorticity and
            ! divergence in spectral space

            complex, dimension(mx,nx), intent(in) :: vorm, divm
            real, dimension(ix,il), intent(out) :: um, vm
            complex, dimension(mx,nx) :: ucosm, vcosm
            complex, dimension(mx,nx) :: zc, zp

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
        end subroutine uvspec
        !*******************************************************************
        subroutine grid(vorm,vorg,kcos)
            ! Transform from spectral to gridpoint space

            real, intent(out) :: vorg(ix,il)
            complex, intent(in) :: vorm(mx,nx)
            integer, intent(in) :: kcos
            complex :: varm(mx,il)

            call gridy(vorm,varm)
            call gridx(varm,vorg,kcos)
        end subroutine grid
        !*********************************************************************
        subroutine spec(vorg,vorm)
            ! Transform from gridpoint to spectral space

            real, intent(in) :: vorg(ix,il)
            complex, intent(out) :: vorm(mx,nx)
            complex :: varm(mx,il)

            call specx(vorg,varm)
            call specy(varm,vorm)
        end subroutine spec
        !*********************************************************************
        subroutine vdspec(ug,vg,vorm,divm,kcos)
            ! Calculate vorticity and divergence in spectral space from u and v
            ! in gridpoint space

            real, intent(in) :: ug(ix,il), vg(ix,il)
            complex, intent(out) :: vorm(mx,nx), divm(mx,nx)
            integer, intent(in) :: kcos
            integer :: i, j
            real :: ug1(ix,il), vg1(ix,il)
            complex :: um(mx,il), vm(mx,il), dumc1(mx,nx), dumc2(mx,nx)

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
        end subroutine vdspec
end module spectral
!*********************************************************************
subroutine gridy(v,varm)
    use mod_atparam
    use spectral, only: cpol, nsh2

    implicit none

    real, intent(in) :: v(mx2,nx)
    real, intent(inout) :: varm(mx2,il)
    real :: vm1(mx2),vm2(mx2)

    integer :: j, j1, m, n

    do j=1,iy
        j1=il+1-j

        do m=1,mx2
            vm1(m)=0.
            vm2(m)=0.
        end do

        do n=1,nx,2
            do m=1,nsh2(n)
                vm1(m)=vm1(m)+v(m,n)*cpol(m,n,j)
            end do
        end do

        do n=2,nx,2
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

    implicit none

    real, intent(in) :: varm(mx2,il)
    real, intent(inout) :: vorm(mx2,nx)
    real :: svarm(mx2,iy), dvarm(mx2,iy)

    integer :: j, j1, m, n

    vorm = 0.0

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
            do m=1,nsh2(n)
                vorm(m,n) = vorm(m,n)+cpol(m,n,j)*svarm(m,j)
            end do
        end do

        do n=2,ntrun1,2
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

    implicit none

    complex, intent(inout) :: vor(mx,nx)

    vor = vor * trfilt
end
