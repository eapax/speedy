subroutine rffti(n, wsave)
!    ****************************************************************
!
!    subroutine rffti(n,wsave)
!
!    ****************************************************************
!
!    subroutine rffti initializes the array wsave which is used in
!    both rfftf and rfftb. the prime factorization of n together with
!    a tabulation of the trigonometric functions are computed and
!    stored in wsave.
!
!    input parameter
!
!    n       the length of the sequence to be transformed.
!
!    output parameter
!
!    wsave   a work array which must be dimensioned at least 2*n+15.
!            the same work array can be used for both rfftf and rfftb
!            as long as n remains unchanged. different wsave arrays
!            are required for different values of n. the contents of
!            wsave must not be changed between calls of rfftf or rfftb.
    use rp_emulator

    implicit none

    integer, intent(in) :: n
    type(rpe_var), intent(inout) :: wsave(*)

    !***first executable statement  rffti
    if (n==1) return
    call rffti1(n,wsave(n+1),wsave(2*n+1))
    return
end subroutine rffti

subroutine rfftb(n, r, wsave)
!    ******************************************************************
!
!    subroutine rfftb(n,r,wsave)
!
!    ******************************************************************
!
!    subroutine rfftb computes the real perodic sequence from its
!    fourier coefficients (fourier synthesis). the transform is defined
!    below at output parameter r.
!
!    input parameters
!
!    n       the length of the array r to be transformed.  the method
!            is most efficient when n is a product of small primes.
!            n may change so long as different work arrays are provided
!
!    r       a real array of length n which contains the sequence
!            to be transformed
!
!    wsave   a work array which must be dimensioned at least 2*n+15.
!            in the program that calls rfftb. the wsave array must be
!            initialized by calling subroutine rffti(n,wsave) and a
!            different wsave array must be used for each different
!            value of n. this initialization does not have to be
!            repeated so long as n remains unchanged thus subsequent
!            transforms can be obtained faster than the first.
!            the same wsave array can be used by rfftf and rfftb.
!
!
!    output parameters
!
!    r       for n even and for i = 1,...,n
!
!                 r(i) = r(1)+(-1)**(i-1)*r(n)
!
!                      plus the sum from k=2 to k=n/2 of
!
!                       2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                      -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!            for n odd and for i = 1,...,n
!
!                 r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
!
!                      2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
!
!                     -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
!
!     *****  note
!                 this transform is unnormalized since a call of rfftf
!                 followed by a call of rfftb will multiply the input
!                 sequence by n.
!
!    wsave   contains results which must not be destroyed between
!            calls of rfftb or rfftf.

    use rp_emulator

    implicit none

    integer, intent(in) :: n
    type(rpe_var), intent(inout) :: r(*), wsave(*)

    !***first executable statement  rfftb
    if (n==1) return
    call rfftb1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
end subroutine rfftb

subroutine rfftf(n, r, wsave)
!    ******************************************************************
!
!    subroutine rfftf(n,r,wsave)
!
!    ******************************************************************
!
!    subroutine rfftf computes the fourier coefficients of a real
!    perodic sequence (fourier analysis). the transform is defined
!    below at output parameter r.
!
!    input parameters
!
!    n       the length of the array r to be transformed.  the method
!            is most efficient when n is a product of small primes.
!            n may change so long as different work arrays are provided
!
!    r       a real array of length n which contains the sequence
!            to be transformed
!
!    wsave   a work array which must be dimensioned at least 2*n+15.
!            in the program that calls rfftf. the wsave array must be
!            initialized by calling subroutine rffti(n,wsave) and a
!            different wsave array must be used for each different
!            value of n. this initialization does not have to be
!            repeated so long as n remains unchanged thus subsequent
!            transforms can be obtained faster than the first.
!            the same wsave array can be used by rfftf and rfftb.
!
!
!    output parameters
!
!    r       r(1) = the sum from i=1 to i=n of r(i)
!
!            if n is even set l =n/2   , if n is odd set l = (n+1)/2
!
!              then for k = 2,...,l
!
!                 r(2*k-2) = the sum from i = 1 to i = n of
!
!                      r(i)*cos((k-1)*(i-1)*2*pi/n)
!
!                 r(2*k-1) = the sum from i = 1 to i = n of
!
!                     -r(i)*sin((k-1)*(i-1)*2*pi/n)
!
!            if n is even
!
!                 r(n) = the sum from i = 1 to i = n of
!
!                      (-1)**(i-1)*r(i)
!
!     *****  note
!                 this transform is unnormalized since a call of rfftf
!                 followed by a call of rfftb will multiply the input
!                 sequence by n.
!
!    wsave   contains results which must not be destroyed between
!            calls of rfftf or rfftb.

    use rp_emulator

    implicit none

    integer, intent(in) :: n
    type(rpe_var), intent(inout) :: r(*), wsave(*)

    !***first executable statement  rfftf
    if (n==1) return
    call rfftf1 (n,r,wsave,wsave(n+1),wsave(2*n+1))
    return
end subroutine rfftf

subroutine rffti1(n, wa, ifac)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: n
    type(rpe_var), intent(inout) :: wa(*)
    integer, intent(inout) :: ifac(*)
    integer, save :: ntryh(4) = (/ 4, 2, 3, 5 /)
    integer :: nl, nf, i, j, ib, ido, ii, ip, ipm, is, k1, l1, l2, ld, nfm1,&
        & nq, nr, ntry
    type(rpe_var) :: arg, argh, argld, fi, tpi

    !***first executable statement  rffti1
      nl = n
      nf = 0
      j = 0
  101 j = j+1
      if (j-4) 102,102,103
  102 ntry = ntryh(j)
      go to 104
  103 ntry = ntry+2
  104 nq = nl/ntry
      nr = nl-ntry*nq
      if (nr) 101,105,101
  105 nf = nf+1
      ifac(nf+2) = ntry
      nl = nq
      if (ntry/=2) go to 107
      if (nf==1) go to 107
      do i=2,nf
         ib = nf-i+2
         ifac(ib+2) = ifac(ib+1)
      end do
      ifac(3) = 2
  107 if (nl  /= 1) go to 104
      ifac(1) = n
      ifac(2) = nf
      tpi = rpe_literal(8.0_dp)*atan(rpe_literal(1.0_dp))
      argh = tpi/n
      is = 0
      nfm1 = nf-1
      l1 = 1
      if (nfm1==0) return
      do k1=1,nfm1
         ip = ifac(k1+2)
         ld = 0
         l2 = l1*ip
         ido = n/l2
         ipm = ip-1
         do j=1,ipm
            ld = ld+l1
            i = is
            argld = ld*argh
            fi = 0.0_dp
            do ii=3,ido,2
               i = i+2
               fi = fi+rpe_literal(1.0_dp)
               arg = fi*argld
               wa(i-1) = cos(arg)
               wa(i) = sin(arg)
            end do
            is = is+ido
         end do
         l1 = l2
      end do
end subroutine rffti1

subroutine rfftb1(n, c, ch, wa, ifac)
    use rp_emulator

    implicit none

    integer, intent(in) :: n, ifac(*)
    type(rpe_var), intent(inout) :: ch(*), c(*), wa(*)
    integer :: nf, na, l1, iw, ip, l2, ido, idl1, ix2, ix3, ix4, i, k1

    !***first executable statement  rfftb1
      nf = ifac(2)
      na = 0
      l1 = 1
      iw = 1
      do k1=1,nf
         ip = ifac(k1+2)
         l2 = ip*l1
         ido = n/l2
         idl1 = ido*l1
         if (ip/=4) go to 103
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na/=0) go to 101
         call radb4(ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 102
  101    call radb4(ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
  102    na = 1-na
         go to 115
  103    if (ip  /= 2) go to 106
         if (na/=0) go to 104
         call radb2(ido,l1,c,ch,wa(iw))
         go to 105
  104    call radb2(ido,l1,ch,c,wa(iw))
  105    na = 1-na
         go to 115
  106    if (ip  /= 3) go to 109
         ix2 = iw+ido
         if (na/=0) go to 107
         call radb3(ido,l1,c,ch,wa(iw),wa(ix2))
         go to 108
  107    call radb3(ido,l1,ch,c,wa(iw),wa(ix2))
  108    na = 1-na
         go to 115
  109    if (ip  /= 5) go to 112
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na/=0) go to 110
         call radb5(ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 111
  110    call radb5(ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
  111    na = 1-na
         go to 115
  112    if (na  /= 0) go to 113
         call radbg(ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         go to 114
  113    call radbg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
  114    if (ido==1) na = 1-na
  115    l1 = l2
         iw = iw+(ip-1)*ido
      end do
      if (na==0) return
      do i=1,n
         c(i) = ch(i)
      end do
end subroutine rfftb1

subroutine rfftf1 (n, c, ch, wa, ifac)
    use rp_emulator

    implicit none

    integer, intent(in) :: n, ifac(*)
    type(rpe_var), intent(inout) :: ch(*), wa(*)
    type(rpe_var), intent(inout) :: c(*)
    integer :: nf, na, l2, iw, k1, kh, ip, l1, ido, idl1, ix2, ix3, ix4, i

    !***FIRST EXECUTABLE STATEMENT  RFFTF1
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip/=4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na/=0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip  /= 2) go to 104
         if (na/=0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip  /= 3) go to 106
         ix2 = iw+ido
         if (na/=0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip  /= 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na/=0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido==1) na = 1-na
         if (na/=0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
      end do
      if (na==1) return
      do i=1,n
         c(i) = ch(i)
      end do
end subroutine rfftf1

subroutine radb2 (ido, l1, cc, ch, wa1)
    use rp_emulator

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,2,*)
    type(rpe_var), intent(inout) :: ch(ido,l1,2), wa1(*)
    integer :: k, idp2, ic, i
    type(rpe_var) :: tr2, ti2

    !***FIRST EXECUTABLE STATEMENT  RADB2
      do k=1,l1
         ch(1,k,1) = cc(1,1,k)+cc(ido,2,k)
         ch(1,k,2) = cc(1,1,k)-cc(ido,2,k)
      end do
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2<l1) go to 108
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
         end do
      end do
      go to 111
  108 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            ch(i-1,k,1) = cc(i-1,1,k)+cc(ic-1,2,k)
            tr2 = cc(i-1,1,k)-cc(ic-1,2,k)
            ch(i,k,1) = cc(i,1,k)-cc(ic,2,k)
            ti2 = cc(i,1,k)+cc(ic,2,k)
            ch(i-1,k,2) = wa1(i-2)*tr2-wa1(i-1)*ti2
            ch(i,k,2) = wa1(i-2)*ti2+wa1(i-1)*tr2
         end do
      end do
  111 if (mod(ido,2)==1) return
  105 do k=1,l1
         ch(ido,k,1) = cc(ido,1,k)+cc(ido,1,k)
         ch(ido,k,2) = -(cc(1,2,k)+cc(1,2,k))
      end do
  107 return
end subroutine radb2

subroutine radb3(ido, l1, cc, ch, wa1, wa2)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,3,*), wa1(*), wa2(*)
    type(rpe_var), intent(inout) :: ch(ido,l1,3)
    integer :: idp2, ic, i, k
    type(rpe_var) :: taur, taui, tr2, cr2, ci3, ti2, ci2, cr3, dr2, dr3, di2, di3

    !***FIRST EXECUTABLE STATEMENT  RADB3
      taur = -0.5_dp
      taui = rpe_literal(0.5_dp)*sqrt(rpe_literal(3.0_dp))
      do k=1,l1
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         cr2 = cc(1,1,k)+taur*tr2
         ch(1,k,1) = cc(1,1,k)+tr2
         ci3 = taui*(cc(1,3,k)+cc(1,3,k))
         ch(1,k,2) = cr2-ci3
         ch(1,k,3) = cr2+ci3
      end do
      if (ido==1) return
      idp2 = ido+2
      if((ido-1)/2<l1) go to 104
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
         end do
      end do
      return
  104 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            cr2 = cc(i-1,1,k)+taur*tr2
            ch(i-1,k,1) = cc(i-1,1,k)+tr2
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ci2 = cc(i,1,k)+taur*ti2
            ch(i,k,1) = cc(i,1,k)+ti2
            cr3 = taui*(cc(i-1,3,k)-cc(ic-1,2,k))
            ci3 = taui*(cc(i,3,k)+cc(ic,2,k))
            dr2 = cr2-ci3
            dr3 = cr2+ci3
            di2 = ci2+cr3
            di3 = ci2-cr3
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
         end do
      end do
end subroutine radb3

subroutine radb4(ido, l1, cc, ch, wa1, wa2, wa3)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,4,*), wa1(*), wa2(*), wa3(*)
    type(rpe_var), intent(inout) :: ch(ido,l1,4)
    type(rpe_var) :: sqrt2, tr1, tr2, tr3, tr4, ti1, ti2, ti3, ti4, cr3, ci3, cr2, cr4,&
        & ci2, ci4
    integer :: i, k, idp2, ic

    !***First executable statement  radb4
      sqrt2 = sqrt(rpe_literal(2.0_dp))
      do k=1,l1
         tr1 = cc(1,1,k)-cc(ido,4,k)
         tr2 = cc(1,1,k)+cc(ido,4,k)
         tr3 = cc(ido,2,k)+cc(ido,2,k)
         tr4 = cc(1,3,k)+cc(1,3,k)
         ch(1,k,1) = tr2+tr3
         ch(1,k,2) = tr1-tr4
         ch(1,k,3) = tr2-tr3
         ch(1,k,4) = tr1+tr4
      end do
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2<l1) go to 108
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
         end do
      end do
      go to 111
  108 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            ti1 = cc(i,1,k)+cc(ic,4,k)
            ti2 = cc(i,1,k)-cc(ic,4,k)
            ti3 = cc(i,3,k)-cc(ic,2,k)
            tr4 = cc(i,3,k)+cc(ic,2,k)
            tr1 = cc(i-1,1,k)-cc(ic-1,4,k)
            tr2 = cc(i-1,1,k)+cc(ic-1,4,k)
            ti4 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr3 = cc(i-1,3,k)+cc(ic-1,2,k)
            ch(i-1,k,1) = tr2+tr3
            cr3 = tr2-tr3
            ch(i,k,1) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1-tr4
            cr4 = tr1+tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(i-1,k,2) = wa1(i-2)*cr2-wa1(i-1)*ci2
            ch(i,k,2) = wa1(i-2)*ci2+wa1(i-1)*cr2
            ch(i-1,k,3) = wa2(i-2)*cr3-wa2(i-1)*ci3
            ch(i,k,3) = wa2(i-2)*ci3+wa2(i-1)*cr3
            ch(i-1,k,4) = wa3(i-2)*cr4-wa3(i-1)*ci4
            ch(i,k,4) = wa3(i-2)*ci4+wa3(i-1)*cr4
         end do
      end do
  111 if (mod(ido,2)==1) return
  105 do k=1,l1
         ti1 = cc(1,2,k)+cc(1,4,k)
         ti2 = cc(1,4,k)-cc(1,2,k)
         tr1 = cc(ido,1,k)-cc(ido,3,k)
         tr2 = cc(ido,1,k)+cc(ido,3,k)
         ch(ido,k,1) = tr2+tr2
         ch(ido,k,2) = sqrt2*(tr1-ti1)
         ch(ido,k,3) = ti2+ti2
         ch(ido,k,4) = -sqrt2*(tr1+ti1)
      end do
  107 return
end subroutine radb4

subroutine radb5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,5,*), wa1(*), wa2(*), wa3(*), wa4(*)
    type(rpe_var), intent(inout) :: ch(ido,l1,5)
    type(rpe_var) :: pi, tr11, ti11, tr12, ti12, ti5, ti4, tr2, tr3, cr2, cr3, ci5, ci4,&
        & ti2, ti3, tr5, tr4, ci2, ci3, cr5, cr4, dr3, dr4, di3, di4, dr5,&
        & dr2, di5, di2
    integer :: i, k, ic, idp2

    !***First executable statement  radb5
      pi = rpe_literal(4.0_dp)*atan(rpe_literal(1.0_dp))
      tr11 = sin(rpe_literal(0.1_dp)*pi)
      ti11 = sin(rpe_literal(0.4_dp)*pi)
      tr12 = -sin(rpe_literal(0.3_dp)*pi)
      ti12 = sin(rpe_literal(0.2_dp)*pi)
      do k=1,l1
         ti5 = cc(1,3,k)+cc(1,3,k)
         ti4 = cc(1,5,k)+cc(1,5,k)
         tr2 = cc(ido,2,k)+cc(ido,2,k)
         tr3 = cc(ido,4,k)+cc(ido,4,k)
         ch(1,k,1) = cc(1,1,k)+tr2+tr3
         cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
         cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
         ci5 = ti11*ti5+ti12*ti4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2) = cr2-ci5
         ch(1,k,3) = cr3-ci4
         ch(1,k,4) = cr3+ci4
         ch(1,k,5) = cr2+ci5
      end do
      if (ido==1) return
      idp2 = ido+2
      if((ido-1)/2<l1) go to 104
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
         end do
      end do
      return
  104 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            ti5 = cc(i,3,k)+cc(ic,2,k)
            ti2 = cc(i,3,k)-cc(ic,2,k)
            ti4 = cc(i,5,k)+cc(ic,4,k)
            ti3 = cc(i,5,k)-cc(ic,4,k)
            tr5 = cc(i-1,3,k)-cc(ic-1,2,k)
            tr2 = cc(i-1,3,k)+cc(ic-1,2,k)
            tr4 = cc(i-1,5,k)-cc(ic-1,4,k)
            tr3 = cc(i-1,5,k)+cc(ic-1,4,k)
            ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
            ch(i,k,1) = cc(i,1,k)+ti2+ti3
            cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
            ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
            cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
            ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(i-1,k,2) = wa1(i-2)*dr2-wa1(i-1)*di2
            ch(i,k,2) = wa1(i-2)*di2+wa1(i-1)*dr2
            ch(i-1,k,3) = wa2(i-2)*dr3-wa2(i-1)*di3
            ch(i,k,3) = wa2(i-2)*di3+wa2(i-1)*dr3
            ch(i-1,k,4) = wa3(i-2)*dr4-wa3(i-1)*di4
            ch(i,k,4) = wa3(i-2)*di4+wa3(i-1)*dr4
            ch(i-1,k,5) = wa4(i-2)*dr5-wa4(i-1)*di5
            ch(i,k,5) = wa4(i-2)*di5+wa4(i-1)*dr5
         end do
      end do
      return
end subroutine radb5

subroutine radbg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, ip, l1, idl1
    type(rpe_var), intent(in) :: cc(ido,ip,*), wa(*)
    type(rpe_var), intent(inout) :: ch(ido,l1,*), c1(ido,l1,*), c2(idl1,*),&
          & ch2(idl1,*)
    type(rpe_var) :: tpi, arg, dcp, dsp, ar1, ai1, ar1h, ds2, dc2, ar2, ai2, ar2h
    integer :: idp2, nbd, ipp2, ipph, i, j, k, jc, j2, is, idij, ic, ik,&
        & l, lc

    !***First executable statement  radbg
      tpi = rpe_literal(8.0_dp)*atan(rpe_literal(1.0_dp))
      arg = tpi/ip
      dcp = cos(arg)
      dsp = sin(arg)
      idp2 = ido+2
      nbd = (ido-1)/2
      ipp2 = ip+2
      ipph = (ip+1)/2
      if (ido<l1) go to 103
      do k=1,l1
         do i=1,ido
            ch(i,k,1) = cc(i,1,k)
         end do
      end do
      go to 106
  103 do i=1,ido
         do k=1,l1
            ch(i,k,1) = cc(i,1,k)
         end do
      end do
  106 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            ch(1,k,j) = cc(ido,j2-2,k)+cc(ido,j2-2,k)
            ch(1,k,jc) = cc(1,j2-1,k)+cc(1,j2-1,k)
         end do
      end do
      if (ido==1) go to 116
      if (nbd<l1) go to 112
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
            end do
         end do
      end do
      go to 116
  112 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            ic = idp2-i
            do k=1,l1
               ch(i-1,k,j) = cc(i-1,2*j-1,k)+cc(ic-1,2*j-2,k)
               ch(i-1,k,jc) = cc(i-1,2*j-1,k)-cc(ic-1,2*j-2,k)
               ch(i,k,j) = cc(i,2*j-1,k)-cc(ic,2*j-2,k)
               ch(i,k,jc) = cc(i,2*j-1,k)+cc(ic,2*j-2,k)
            end do
         end do
       end do
  116 ar1 = 1.0_dp
      ai1 = 0.0_dp
      do l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do ik=1,idl1
            c2(ik,l) = ch2(ik,1)+ar1*ch2(ik,2)
            c2(ik,lc) = ai1*ch2(ik,ip)
         end do
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do ik=1,idl1
               c2(ik,l) = c2(ik,l)+ar2*ch2(ik,j)
               c2(ik,lc) = c2(ik,lc)+ai2*ch2(ik,jc)
            end do
         end do
      end do
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
         end do
      end do
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            ch(1,k,j) = c1(1,k,j)-c1(1,k,jc)
            ch(1,k,jc) = c1(1,k,j)+c1(1,k,jc)
         end do
      end do
      if (ido==1) go to 132
      if (nbd<l1) go to 128
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
            end do
         end do
      end do
      go to 132
  128 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            do k=1,l1
               ch(i-1,k,j) = c1(i-1,k,j)-c1(i,k,jc)
               ch(i-1,k,jc) = c1(i-1,k,j)+c1(i,k,jc)
               ch(i,k,j) = c1(i,k,j)+c1(i-1,k,jc)
               ch(i,k,jc) = c1(i,k,j)-c1(i-1,k,jc)
            end do
         end do
      end do
  132 continue
      if (ido==1) return
      do ik=1,idl1
         c2(ik,1) = ch2(ik,1)
      end do
      do j=2,ip
         do k=1,l1
            c1(1,k,j) = ch(1,k,j)
         end do
      end do
      if (nbd>l1) go to 139
      is = -ido
      do j=2,ip
         is = is+ido
         idij = is
         do i=3,ido,2
            idij = idij+2
            do k=1,l1
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
            end do
         end do
      end do
      go to 143
  139 is = -ido
      do j=2,ip
         is = is+ido
         do k=1,l1
            idij = is
            do i=3,ido,2
               idij = idij+2
               c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
               c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
            end do
         end do
      end do
  143 return
end subroutine radbg

subroutine radf2(ido, l1, cc, ch, wa1)
    use rp_emulator

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,l1,2), wa1(*)
    type(rpe_var), intent(inout) :: ch(ido,2,*)
    type(rpe_var) :: tr2, ti2
    integer :: i, k, idp2, ic

    !***First executable statement  radf2
      do k=1,l1
         ch(1,1,k) = cc(1,k,1)+cc(1,k,2)
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,2)
      end do
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2<l1) go to 108
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
         end do
      end do
      go to 111
  108 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            tr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ti2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            ch(i,1,k) = cc(i,k,1)+ti2
            ch(ic,2,k) = ti2-cc(i,k,1)
            ch(i-1,1,k) = cc(i-1,k,1)+tr2
            ch(ic-1,2,k) = cc(i-1,k,1)-tr2
         end do
      end do
  111 if (mod(ido,2)==1) return
  105 do k=1,l1
         ch(1,2,k) = -cc(ido,k,2)
         ch(ido,1,k) = cc(ido,k,1)
      end do
  107 return
end subroutine radf2

subroutine radf3(ido, l1, cc, ch, wa1, wa2)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,l1,3), wa1(*), wa2(*)
    type(rpe_var), intent(inout) :: ch(ido,3,*)
    type(rpe_var) :: taur, taui, cr2, dr2, di2, dr3, di3, ci2, tr2, ti2, tr3, ti3
    integer :: i, k, idp2, ic

    !***First executable statement  radf3
      taur = -0.5_dp
      taui = rpe_literal(0.5_dp)*sqrt(rpe_literal(3.0_dp))
      do k=1,l1
         cr2 = cc(1,k,2)+cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2
         ch(1,3,k) = taui*(cc(1,k,3)-cc(1,k,2))
         ch(ido,2,k) = cc(1,k,1)+taur*cr2
      end do
      if (ido==1) return
      idp2 = ido+2
      if((ido-1)/2<l1) go to 104
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
         end do
      end do
      return
  104 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr2 = dr2+dr3
            ci2 = di2+di3
            ch(i-1,1,k) = cc(i-1,k,1)+cr2
            ch(i,1,k) = cc(i,k,1)+ci2
            tr2 = cc(i-1,k,1)+taur*cr2
            ti2 = cc(i,k,1)+taur*ci2
            tr3 = taui*(di2-di3)
            ti3 = taui*(dr3-dr2)
            ch(i-1,3,k) = tr2+tr3
            ch(ic-1,2,k) = tr2-tr3
            ch(i,3,k) = ti2+ti3
            ch(ic,2,k) = ti3-ti2
         end do
      end do
end subroutine radf3

subroutine radf4(ido, l1, cc, ch, wa1, wa2, wa3)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,l1,4), wa1(*), wa2(*), wa3(*)
    type(rpe_var), intent(inout) :: ch(ido,4,*)
    type(rpe_var) :: hsqt2, tr1, tr2, tr3, tr4, cr2, ci2, cr3, ci3, cr4, ci4, ti1, ti2,&
        & ti3, ti4
    integer :: i, k, ic, idp2

    !***First executable statement  radf4
      hsqt2 = rpe_literal(0.5_dp)*sqrt(rpe_literal(2.0_dp))
      do k=1,l1
         tr1 = cc(1,k,2)+cc(1,k,4)
         tr2 = cc(1,k,1)+cc(1,k,3)
         ch(1,1,k) = tr1+tr2
         ch(ido,4,k) = tr2-tr1
         ch(ido,2,k) = cc(1,k,1)-cc(1,k,3)
         ch(1,3,k) = cc(1,k,4)-cc(1,k,2)
      end do
      if (ido-2) 107,105,102
  102 idp2 = ido+2
      if((ido-1)/2<l1) go to 111
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
         end do
      end do
      go to 110
  111 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            cr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            ci2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            cr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            ci3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            cr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            ci4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            tr1 = cr2+cr4
            tr4 = cr4-cr2
            ti1 = ci2+ci4
            ti4 = ci2-ci4
            ti2 = cc(i,k,1)+ci3
            ti3 = cc(i,k,1)-ci3
            tr2 = cc(i-1,k,1)+cr3
            tr3 = cc(i-1,k,1)-cr3
            ch(i-1,1,k) = tr1+tr2
            ch(ic-1,4,k) = tr2-tr1
            ch(i,1,k) = ti1+ti2
            ch(ic,4,k) = ti1-ti2
            ch(i-1,3,k) = ti4+tr3
            ch(ic-1,2,k) = tr3-ti4
            ch(i,3,k) = tr4+ti3
            ch(ic,2,k) = tr4-ti3
         end do
      end do
  110 if (mod(ido,2)==1) return
  105 do k=1,l1
         ti1 = -hsqt2*(cc(ido,k,2)+cc(ido,k,4))
         tr1 = hsqt2*(cc(ido,k,2)-cc(ido,k,4))
         ch(ido,1,k) = tr1+cc(ido,k,1)
         ch(ido,3,k) = cc(ido,k,1)-tr1
         ch(1,2,k) = ti1-cc(ido,k,3)
         ch(1,4,k) = ti1+cc(ido,k,3)
      end do
  107 return
end subroutine radf4

subroutine radf5(ido, l1, cc, ch, wa1, wa2, wa3, wa4)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, l1
    type(rpe_var), intent(in) :: cc(ido,l1,5), wa1(*), wa2(*), wa3(*), wa4(*)
    type(rpe_var), intent(inout) :: ch(ido,5,*)
    type(rpe_var) :: pi, tr11, ti11, tr12, ti12, cr2, cr3, cr4, cr5, ci2, ci3, ci4, ci5,&
        & dr2, dr3, dr4, dr5, di2, di3, di4, di5, tr2, tr3, tr4, tr5, ti2, ti3,&
        & ti4, ti5
    integer :: i, k, idp2, ic

    !***First executable statement  radf5
      pi = rpe_literal(4.0_dp)*atan(rpe_literal(1.0_dp))
      tr11 = sin(rpe_literal(0.1_dp)*pi)
      ti11 = sin(rpe_literal(0.4_dp)*pi)
      tr12 = -sin(rpe_literal(0.3_dp)*pi)
      ti12 = sin(rpe_literal(0.2_dp)*pi)
      do k=1,l1
         cr2 = cc(1,k,5)+cc(1,k,2)
         ci5 = cc(1,k,5)-cc(1,k,2)
         cr3 = cc(1,k,4)+cc(1,k,3)
         ci4 = cc(1,k,4)-cc(1,k,3)
         ch(1,1,k) = cc(1,k,1)+cr2+cr3
         ch(ido,2,k) = cc(1,k,1)+tr11*cr2+tr12*cr3
         ch(1,3,k) = ti11*ci5+ti12*ci4
         ch(ido,4,k) = cc(1,k,1)+tr12*cr2+tr11*cr3
         ch(1,5,k) = ti12*ci5-ti11*ci4
      end do
      if (ido==1) return
      idp2 = ido+2
      if((ido-1)/2<l1) go to 104
      do k=1,l1
         do i=3,ido,2
            ic = idp2-i
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
         end do
      end do
      return
  104 do i=3,ido,2
         ic = idp2-i
         do k=1,l1
            dr2 = wa1(i-2)*cc(i-1,k,2)+wa1(i-1)*cc(i,k,2)
            di2 = wa1(i-2)*cc(i,k,2)-wa1(i-1)*cc(i-1,k,2)
            dr3 = wa2(i-2)*cc(i-1,k,3)+wa2(i-1)*cc(i,k,3)
            di3 = wa2(i-2)*cc(i,k,3)-wa2(i-1)*cc(i-1,k,3)
            dr4 = wa3(i-2)*cc(i-1,k,4)+wa3(i-1)*cc(i,k,4)
            di4 = wa3(i-2)*cc(i,k,4)-wa3(i-1)*cc(i-1,k,4)
            dr5 = wa4(i-2)*cc(i-1,k,5)+wa4(i-1)*cc(i,k,5)
            di5 = wa4(i-2)*cc(i,k,5)-wa4(i-1)*cc(i-1,k,5)
            cr2 = dr2+dr5
            ci5 = dr5-dr2
            cr5 = di2-di5
            ci2 = di2+di5
            cr3 = dr3+dr4
            ci4 = dr4-dr3
            cr4 = di3-di4
            ci3 = di3+di4
            ch(i-1,1,k) = cc(i-1,k,1)+cr2+cr3
            ch(i,1,k) = cc(i,k,1)+ci2+ci3
            tr2 = cc(i-1,k,1)+tr11*cr2+tr12*cr3
            ti2 = cc(i,k,1)+tr11*ci2+tr12*ci3
            tr3 = cc(i-1,k,1)+tr12*cr2+tr11*cr3
            ti3 = cc(i,k,1)+tr12*ci2+tr11*ci3
            tr5 = ti11*cr5+ti12*cr4
            ti5 = ti11*ci5+ti12*ci4
            tr4 = ti12*cr5-ti11*cr4
            ti4 = ti12*ci5-ti11*ci4
            ch(i-1,3,k) = tr2+tr5
            ch(ic-1,2,k) = tr2-tr5
            ch(i,3,k) = ti2+ti5
            ch(ic,2,k) = ti5-ti2
            ch(i-1,5,k) = tr3+tr4
            ch(ic-1,4,k) = tr3-tr4
            ch(i,5,k) = ti3+ti4
            ch(ic,4,k) = ti4-ti3
         end do
      end do
end subroutine radf5

subroutine radfg(ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa)
    use rp_emulator
    use mod_prec

    implicit none

    integer, intent(in) :: ido, ip, l1, idl1
    type(rpe_var), intent(in) :: wa(*)
    type(rpe_var), intent(inout) :: ch(ido,l1,*), cc(ido,ip,*), c1(ido,l1,*),&
        & c2(idl1,*), ch2(idl1,*)
    type(rpe_var) :: tpi, arg, dcp, dsp, ar1h, ar2h, ai1, ai2, ar1, ar2, dc2, ds2
    integer :: ipph, ipp2, idp2, nbd, is, ik, j, j2, jc, i, ic, idij, k,&
        & l, lc

    !***First executable statement  radfg
      tpi = rpe_literal(8.0_dp)*atan(rpe_literal(1.0_dp))
      arg = tpi/ip
      dcp = cos(arg)
      dsp = sin(arg)
      ipph = (ip+1)/2
      ipp2 = ip+2
      idp2 = ido+2
      nbd = (ido-1)/2
      if (ido==1) go to 119
      do ik=1,idl1
         ch2(ik,1) = c2(ik,1)
      end do
      do j=2,ip
         do k=1,l1
            ch(1,k,j) = c1(1,k,j)
         end do
      end do
      if (nbd>l1) go to 107
      is = -ido
      do j=2,ip
         is = is+ido
         idij = is
         do i=3,ido,2
            idij = idij+2
            do k=1,l1
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
            end do
         end do
      end do
      go to 111
  107 is = -ido
      do j=2,ip
         is = is+ido
         do k=1,l1
            idij = is
            do i=3,ido,2
               idij = idij+2
               ch(i-1,k,j) = wa(idij-1)*c1(i-1,k,j)+wa(idij)*c1(i,k,j)
               ch(i,k,j) = wa(idij-1)*c1(i,k,j)-wa(idij)*c1(i-1,k,j)
            end do
         end do
      end do
  111 if (nbd < l1) go to 115
      do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            do i=3,ido,2
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
            end do
         end do
      end do
      go to 121
  115 do j=2,ipph
         jc = ipp2-j
         do i=3,ido,2
            do k=1,l1
               c1(i-1,k,j) = ch(i-1,k,j)+ch(i-1,k,jc)
               c1(i-1,k,jc) = ch(i,k,j)-ch(i,k,jc)
               c1(i,k,j) = ch(i,k,j)+ch(i,k,jc)
               c1(i,k,jc) = ch(i-1,k,jc)-ch(i-1,k,j)
            end do
         end do
      end do
      go to 121
  119 do ik=1,idl1
         c2(ik,1) = ch2(ik,1)
      end do
  121 do j=2,ipph
         jc = ipp2-j
         do k=1,l1
            c1(1,k,j) = ch(1,k,j)+ch(1,k,jc)
            c1(1,k,jc) = ch(1,k,jc)-ch(1,k,j)
         end do
      end do
!
      ar1 = 1.0_dp
      ai1 = 0.0_dp
      do l=2,ipph
         lc = ipp2-l
         ar1h = dcp*ar1-dsp*ai1
         ai1 = dcp*ai1+dsp*ar1
         ar1 = ar1h
         do ik=1,idl1
            ch2(ik,l) = c2(ik,1)+ar1*c2(ik,2)
            ch2(ik,lc) = ai1*c2(ik,ip)
         end do
         dc2 = ar1
         ds2 = ai1
         ar2 = ar1
         ai2 = ai1
         do j=3,ipph
            jc = ipp2-j
            ar2h = dc2*ar2-ds2*ai2
            ai2 = dc2*ai2+ds2*ar2
            ar2 = ar2h
            do ik=1,idl1
               ch2(ik,l) = ch2(ik,l)+ar2*c2(ik,j)
               ch2(ik,lc) = ch2(ik,lc)+ai2*c2(ik,jc)
            end do
         end do
      end do
      do j=2,ipph
         do ik=1,idl1
            ch2(ik,1) = ch2(ik,1)+c2(ik,j)
         end do
      end do
!
      if (ido<l1) go to 132
      do k=1,l1
         do i=1,ido
            cc(i,1,k) = ch(i,k,1)
         end do
      end do
      go to 135
  132 do i=1,ido
         do k=1,l1
            cc(i,1,k) = ch(i,k,1)
         end do
      end do
  135 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            cc(ido,j2-2,k) = ch(1,k,j)
            cc(1,j2-1,k) = ch(1,k,jc)
         end do
      end do
      if (ido==1) return
      if (nbd<l1) go to 141
      do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do k=1,l1
            do i=3,ido,2
               ic = idp2-i
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
            end do
         end do
      end do
      return
  141 do j=2,ipph
         jc = ipp2-j
         j2 = j+j
         do i=3,ido,2
            ic = idp2-i
            do k=1,l1
               cc(i-1,j2-1,k) = ch(i-1,k,j)+ch(i-1,k,jc)
               cc(ic-1,j2-2,k) = ch(i-1,k,j)-ch(i-1,k,jc)
               cc(i,j2-1,k) = ch(i,k,j)+ch(i,k,jc)
               cc(ic,j2-2,k) = ch(i,k,jc)-ch(i,k,j)
            end do
         end do
      end do
      return
end subroutine radfg
