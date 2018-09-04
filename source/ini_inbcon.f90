subroutine inbcon(grav0,radlat)
    !
    ! subroutine inbcon (grav0,radlat)
    !
    ! Purpose : Read topography and climatological boundary conditions
    ! Input :   grav0  = gravity accel.
    !           radlat = grid latitudes in radiants

    use mod_cpl_flags, only: icsea, isstan
    use mod_date, only: isst0
    use mod_atparam
    use mod_surfcon
    use mod_cli_land
    use mod_cli_sea
    use mod_prec, only: sp, dp

    implicit none

    real(dp), intent(in) :: grav0
    real(dp), intent(in) :: radlat(il)

    real(sp) :: r4inp(ix,il)
    real(sp) :: inp(ix,il)
    real(sp) :: veg(ix,il), swl1(ix,il), swl2(ix,il)

    integer :: iitest=1, i, idep2, irec, irecl, it, j, jrec
    real(dp) :: rad2deg, rsw, sdep1, sdep2, swroot, swwil2, thrsh

    ! Set threshold for land-sea mask definition
    ! (ie minimum fraction of either land or sea)

    thrsh = 0.1_dp

    ! 1. Read topographical fields (orography, land-sea mask)
    if (iitest >= 1) print *,' read orography'

    call load_boundary_file(20,inp,0)

    phi0 = grav0*inp

    call truncg (ntrun,phi0,phis0)

    if (iitest >= 1) print *,' read fractional land-sea mask'

    call load_boundary_file(20,inp,1)

    fmask = inp

    ! 2. Initialize land-sfc boundary conditions

    ! 2.1 Fractional and binary land masks
    do j=1,il
        do i=1,ix
            fmask_l(i,j) = fmask(i,j)

            if (fmask_l(i,j).ge.thrsh) then
              bmask_l(i,j) = 1.0_dp
              if (fmask(i,j).gt.(1.0_dp-thrsh)) fmask_l(i,j) = 1.0_dp
            else
              bmask_l(i,j) = 0.0_dp
              fmask_l(i,j) = 0.0_dp
            end if

            fmask1(i,j) = fmask_l(i,j)
        end do
    end do

    ! 2.2 Annual-mean surface albedo
    if (iitest >= 1) print *,' read surface albedo'

    call load_boundary_file(20,inp,2)

    alb0 = inp

    ! 2.3 Land-surface temp.
    if (iitest >= 1) print *,' reading land-surface temp.'

    do it = 1,12
        call load_boundary_file(23,inp,it-1)

        call fillsf(inp,ix,il,0.0_dp)

       stl12(1:ix,1:il,it) = inp
    end do

    if (iitest == 1) print *,' checking land-surface temp.'

    call forchk(bmask_l,stl12,ngp,12,0.0_dp,400.0_dp,273.0_dp)

    ! 2.4 Snow depth
    if (iitest >= 1) print *,' reading snow depth'

    do it = 1,12
        call load_boundary_file(24,inp,it-1)

        snowd12(1:ix,1:il,it) = inp
    end do

    if (iitest >= 1) print *,' checking snow depth'

    CALL FORCHK (bmask_l,snowd12,ngp,12,0.0_dp,20000.0_dp,0.0_dp)

    ! 2.5 Read soil moisture and compute soil water availability
    !     using vegetation fraction

    if (iitest >= 1) print *,' reading soil moisture'

    ! Read vegetation fraction
    call load_boundary_file(20,veg,3)
    call load_boundary_file(20,inp,4)

    ! Combine high and low vegetation fractions
    veg = max(0.0_dp,veg+0.8_dp*inp)

    ! Read soil moisture
    sdep1 = 70.0_dp
    idep2 = 3
    sdep2 = idep2*sdep1

    swwil2= idep2*swwil
    rsw   = 1.0_dp/(swcap+idep2*(swcap-swwil))

    do it = 1,12
        call load_boundary_file(26,swl1,3*it-3)
        call load_boundary_file(26,swl2,3*it-2)

        ! Combine soil water content from two top layers
        do j = 1,il
            do i = 1,ix
                swroot = idep2*swl2(i,j)
                inp(i,j) = min(1.0_dp, rsw*(swl1(i,j) + &
                        veg(i,j)*max(0.0_dp, swroot-swwil2)))
            end do
        end do

        soilw12(1:ix,1:il,it) = inp
    end do

    if (iitest >= 1) print *,' checking soil moisture'

    call forchk(bmask_l,soilw12,ngp,12,0.0_dp,10.0_dp,0.0_dp)

    ! 3. Initialize sea-sfc boundary conditions

    ! 3.1 Fractional and binary sea masks
    do j=1,il
        do i=1,ix
            fmask_s(i,j) = 1.0_dp-fmask(i,j)

            if (fmask_s(i,j).ge.thrsh) then
                bmask_s(i,j) = 1.0_dp
                if (fmask_s(i,j).gt.(1.0_dp-thrsh)) fmask_s(i,j) = 1.0_dp
            else
                bmask_s(i,j) = 0.0_dp
                fmask_s(i,j) = 0.0_dp
            end if
        end do
    end do

    ! Grid latitudes for sea-sfc. variables
    rad2deg = 90.0_dp/asin(1.0_dp)
    deglat_s = rad2deg*radlat

    ! 3.2 SST
    if (iitest >= 1) print *,' reading sst'

    do it = 1,12
        call load_boundary_file(21,inp,it-1)

        call fillsf(inp,ix,il,0.0_dp)

        sst12(1:ix,1:il,it) = inp
    end do

    if (iitest >= 1) print *,' checking sst'

    call forchk(bmask_s,sst12,ngp,12,100.0_dp,400.0_dp,273.0_dp)

    ! 3.3 Sea ice concentration
    if (iitest >= 1) print *,' reading sea ice'

    do it = 1,12
        call load_boundary_file(22,inp,it-1)

        inp = max(inp,0.0_dp)

        sice12(1:ix,1:il,it) = inp
    end do

    if (iitest >= 1) print *,' checking sea ice'

    call forchk(bmask_s,sice12,ngp,12,0.0_dp,1.0_dp,0.0_dp)

    ! 3.4 SST anomalies for initial and prec./following months
    if (isstan > 0) then
        if (iitest >= 1) print *,' reading sst anomalies'

        print *, 'isst0 = ', isst0
        do it=1,3
            if ((isst0 <= 1 .and. it /= 2) .or. isst0 > 1) then
                call load_boundary_file(30,inp,isst0-2+it-1)
            end if

            sstan3(1:ix,1:il,it) = inp
        end do

        if (iitest >= 1) print *,' checking sst anomalies'

        call forchk(bmask_s,sstan3,ngp,3,-50.0_dp,50.0_dp,0.0_dp)
    end if

    ! 4. Climatological fields for the ocean model (TO BE RECODED)
    ! 4.1. Annual-mean heat flux into sea-surface

    hfseacl = 0.0_dp

    if (icsea >= 1) then
        if (iitest >= 1) print *,' reading sfc heat fluxes'

        irecl = 4*ngp
        irec = 0

        open ( unit=31, file='fort.31', status='old',&
            & form='unformatted', access='direct', recl=irecl )

        do it = 1,12
            irec=irec+2
            read (31,rec=irec) r4inp

            do j = 1,il
                do i = 1,ix
                    hfseacl(i,j) = hfseacl(i,j)+r4inp(i,j)
                end do
            end do
        end do

        do j = 1,il
            do i = 1,ix
                if (bmask_s(i,j).gt.0.0_dp) then
                    hfseacl(i,j) = hfseacl(i,j)/(12.0_dp*fmask_s(i,j))
                else
                    hfseacl(i,j) = 0.0_dp
                end if
            end do
        end do

        if (iitest >= 1) print *,' checking sfc heat fluxes'

        call forchk (bmask_s,hfseacl,ngp,1,-1000.0_dp,1000.0_dp,0.0_dp)
    end if

    ! 4.2. Ocean model SST climatology:
    !      defined by adding SST model bias to obs. climatology
    !      (bias may be defined in a different period from climatology)

    if (icsea >= 3) then
        if (iitest >= 1) print *,' reading ocean model SST bias'

        do it = 1,12
            read (32) r4inp

            do j = 1,il
                do i = 1,ix
                    sstom12(i,j,it) = sst12(i,j,it)+r4inp(i,j)
                end do
            end do
        end do

        if (iitest >= 1) print *,' checking ocean model SST'

        call forchk (bmask_s,sstom12,ngp,12,100.0_dp,400.0_dp,273.0_dp)
    end if
end

subroutine forchk (fmask,field,ngp,nf,fmin,fmax,fset)
    ! Aux. routine forchk: Check consistency of sfc fields with land-sea mask
    ! and set undefined values to a constant (to avoid over/underflow)

    use mod_prec, only: dp

    implicit none

    real(dp), intent(in) :: fmask(ngp)
    real(dp), intent(inout) :: field(ngp,nf)
    integer, intent(in) :: ngp, nf
    real(dp), intent(in) :: fmin, fmax, fset

    integer :: jf, jgp, nfault

    do jf = 1,nf
        nfault=0

        do jgp = 1,ngp
            if (fmask(jgp).gt.0.0_dp) then
                if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax) then
                    nfault = nfault+1
                end if
            else
                field(jgp,jf) = fset
            end if
        end do

        print *, ' field: ', jf, '   no. of faulty points:', nfault
    end do

    print *, ' undefined values set to', fset
end

subroutine truncg (itr,fg1,fg2)
    ! subroutine truncg (itr,fg1,fg2)
    ! Purpose : compute a spectrally-filtered grid-point field
    ! Input   : itr : spectral truncation (triangular)
    !         : fg1 : original grid-point field
    ! Output  : fg2 : filtered grid-point field

    USE mod_atparam
    use spectral, only: grid, spec
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: itr

    real(dp), intent(in) :: fg1 (ix,il)
    real(dp), intent(out) :: fg2(ix,il)
    complex(dp) :: fsp(mx,nx), zero
    integer :: n, m, itwn

    print *, 'Filter applied at wavenumber ', itr

    zero = (0.0_dp,0.0_dp)

    call spec (fg1,fsp)

    do n=1,nx
        do m=1,mx
            itwn=isc*(m-1)+n-1
            if (itwn.gt.itr) fsp(m,n)=zero
        end do
    end do

    call grid (fsp,fg2,1)
end

subroutine fillsf(sf,ix,il,fmis)
    ! subroutine fillsf (sf,ix,nlat)
    ! Purpose: replace missing values in surface fields
    ! NB: it is assumed that non-missing values exist near the Equator

    use mod_prec, only: sp, dp

    implicit none

    real(sp), intent(inout) :: sf(ix,il)
    integer, intent(in) :: ix, il
    real(dp), intent(in) :: fmis

    real(dp) :: sf2(0:ix+1)
    integer :: khem, j, j1, j2, j3, i, nmis
    real(dp) :: fmean

    do khem = 1,2
       if (khem == 1) then
            j1 = il/2
            j2 = 1
            j3 = -1
        else
            j1 = j1+1
            j2 = il
            j3 = 1
        end if

        do j=j1,j2,j3
            sf2(1:ix) = sf(1:ix,j)

            nmis = 0
            do i=1,ix
                if (sf(i,j) < fmis) then
                    nmis = nmis+1
                    sf2(i) = 0.0_dp
                end if
            end do

            if (nmis < ix) fmean = sum(sf2(1:ix))/float(ix-nmis)

            do i=1,ix
                if (sf(i,j).lt.fmis) sf2(i) = fmean
            end do

            sf2(0)      = sf2(ix)
            sf2(ix+1) = sf2(1)
            do i=1,ix
                if (sf(i,j).lt.fmis) sf(i,j) = 0.5_dp*(sf2(i-1)+sf2(i+1))
            end do
        end do
    end do
end

subroutine load_boundary_file(iunit,inp,offset)
    ! read field on from unit=iunit

    use mod_atparam
    use mod_prec, only: sp

    implicit none

    integer, intent(in) :: iunit, offset
    real(sp), intent(out) :: inp(ix,il)
    integer :: i

    open(unit=iunit, form='unformatted', access='direct', recl=ix*4, &
            convert='little_endian')
    do i = 1, il
        read(iunit,rec=offset*il+i) inp(:,il+1-i)
    end do

    ! Fix undefined values
    where (inp <= -999) inp = 0.0_sp

    close(unit=iunit)
end
