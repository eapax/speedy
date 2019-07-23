subroutine inbcon(grav0,radlat)
    !
    ! subroutine inbcon (grav0,radlat)
    !
    ! Purpose : Read topography and climatological boundary conditions
    ! Input :   grav0  = gravity accel.
    !           radlat = grid latitudes in radiants

    use netcdf
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
    real(dp) :: vegh(ix,il), vegl(ix,il), veg(ix,il)
    real(dp) :: swl1(ix,il,12), swl2(ix,il,12)

    integer :: iitest=1, i, idep2, irec, irecl, it, j, NCID, VARID
    real(dp) :: rad2deg, rsw, sdep1, sdep2, swwil2, thrsh, swroot

    ! Set threshold for land-sea mask definition
    ! (ie minimum fraction of either land or sea)
    thrsh = 0.1_dp

    ! 1. Read topographical fields (orography, land-sea mask)
    if (iitest>=1) print *,' read orography'

    call check(NF90_OPEN('climatology.nc', NF90_NOWRITE, NCID))
    call NC_extract_variable(NCID, 'orog', 1, phi0)

    phi0 = grav0*phi0

    call truncg (ntrun,phi0,phis0)

    if (iitest>=1) print *,' read fractional land-sea mask'

    call NC_extract_variable(NCID, 'lsm', 1, fmask)

    ! 2. Initialize land-sfc boundary conditions

    ! 2.1 Fractional and binary land masks
    where (fmask > 1.0_dp-thrsh)
        fmask_l = 1.0_dp
        bmask_l = 1.0_dp
    elsewhere (fmask >= thrsh)
        fmask_l = fmask
        bmask_l = 1.0_dp
    elsewhere
        bmask_l = 0.0_dp
        fmask_l = 0.0_dp
    end where

    fmask1 = fmask_l

    ! 2.2 Annual-mean surface albedo
    if (iitest>=1) print *,' read surface albedo'

    call NC_extract_variable(NCID, 'alb', 1, alb0)

    ! 2.3 Land-surface temp.
    if (iitest>=1) print *,' reading land-surface temp.'

    call NC_extract_variable(NCID, 'stl', 12, stl12)

    do it = 1,12
        call fillsf(stl12(:,:,it),ix,il,0.0_dp)
    end do

    if (iitest==1) print *,' checking land-surface temp.'

    call forchk(bmask_l,stl12,ngp,12,0.0_dp,400.0_dp,273.0_dp)

    ! 2.4 Snow depth
    if (iitest>=1) print *,' reading snow depth'
    call NC_extract_variable(NCID, 'snowd', 12, snowd12)

    if (iitest>=1) print *,' checking snow depth'
    call FORCHK (bmask_l,snowd12,ngp,12,0.0_dp,20000.0_dp,0.0_dp)

    ! 2.5 Read soil moisture and compute soil water availability
    !     using vegetation fraction
    if (iitest>=1) print *,' reading soil moisture'

    ! Read vegetation fraction
    call NC_extract_variable(NCID, 'vegh', 1, vegh)
    call NC_extract_variable(NCID, 'vegl', 1, vegl)

    ! Combine high and low vegetation fractions
    veg = max(0.0_dp, vegh+0.8_dp*vegl)

    ! Read soil moisture
    sdep1 = 70.0_dp
    idep2 = 3
    sdep2 = idep2*sdep1

    swwil2= idep2*swwil
    rsw   = 1.0_dp/(swcap+idep2*(swcap-swwil))

    call NC_extract_variable(NCID, 'swl1', 12, swl1)
    call NC_extract_variable(NCID, 'swl2', 12, swl2)
    do it = 1,12
        ! Combine soil water content from two top layers
        do j = 1,il
            do i = 1,ix
                swroot = idep2*swl2(i,j,it)
                inp(i,j) = min(1.0_dp, rsw*(swl1(i,j,it) + &
                        veg(i,j)*max(0.0_dp, swroot-swwil2)))
            end do
        end do

        soilw12(1:ix,1:il,it) = inp
    end do

    if (iitest>=1) print *,' checking soil moisture'
    call forchk(bmask_l,soilw12,ngp,12,0.0_dp,10.0_dp,0.0_dp)

    ! 3. Initialize sea-sfc boundary conditions

    ! 3.1 Fractional and binary sea masks
    where (fmask < thrsh)
        fmask_s = 1.0_dp
        bmask_s = 1.0_dp
    elsewhere (fmask <= 1.0_dp - thrsh)
        fmask_s = 1.0_dp - fmask
        bmask_s = 1.0_dp
    elsewhere
        bmask_s = 0.0_dp
        fmask_s = 0.0_dp
    end where

    ! Grid latitudes for sea-sfc. variables
    rad2deg = 90.0_dp/asin(1.0_dp)
    deglat_s = rad2deg*radlat

    ! 3.2 SST
    if (iitest>=1) print *,' reading sst'
    call NC_extract_variable(NCID, 'sst', 12, sst12)

    do it = 1,12
        call fillsf(sst12(:,:,it),ix,il,0.0_dp)
    end do

    if (iitest>=1) print *,' checking sst'
    call forchk(bmask_s,sst12,ngp,12,100.0_dp,400.0_dp,273.0_dp)

    ! 3.3 Sea ice concentration
    if (iitest>=1) print *,' reading sea ice'
    call NC_extract_variable(NCID, 'icec', 12, sice12)
    sice12 = max(sice12, 0.0_dp)

    if (iitest>=1) print *,' checking sea ice'
    call forchk(bmask_s,sice12,ngp,12,0.0_dp,1.0_dp,0.0_dp)

    call check(NF90_CLOSE(NCID))

    ! 3.4 SST anomalies for initial and prec./following months
    if (isstan>0) then
        if (iitest>=1) print *,' reading sst anomalies'

        print *, 'isst0 = ', isst0
        ! Read in the the sst anomaly for the current month and the month
        ! before and after the current month
        ! isst0 = number of months of the current date from the start of the
        !         sst anomaly file
        call check(NF90_OPEN('anomalies.nc', NF90_NOWRITE, NCID))
        do it=1,3
            call check(NF90_INQ_VARID(NCID, 'ssta', VARID))
            call check(NF90_GET_VAR(NCID, VARID, inp, &
                    start=(/1, 1, 1, isst0-1+it-1/), count=(/ix, il, 1, 1/) ))

            sstan3(1:ix,1:il,it) = inp
        end do
        call check(NF90_CLOSE(NCID))

        if (iitest>=1) print *,' checking sst anomalies'
        call forchk(bmask_s,sstan3,ngp,3,-50.0_dp,50.0_dp,0.0_dp)
    end if

    ! 4. Climatological fields for the ocean model (TO BE RECODED)
    ! 4.1. Annual-mean heat flux into sea-surface

    hfseacl = 0.0_dp

    if (icsea>=1) then
        if (iitest>=1) print *,' reading sfc heat fluxes'

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

        where (bmask_s>0.0_dp)
            hfseacl = hfseacl/(12.0_dp*fmask_s)
        elsewhere
            hfseacl = 0.0_dp
        end where

        if (iitest>=1) print *,' checking sfc heat fluxes'
        call forchk (bmask_s,hfseacl,ngp,1,-1000.0_dp,1000.0_dp,0.0_dp)
    end if

    ! 4.2. Ocean model SST climatology:
    !      defined by adding SST model bias to obs. climatology
    !      (bias may be defined in a different period from climatology)

    if (icsea>=3) then
        if (iitest>=1) print *,' reading ocean model SST bias'

        do it = 1,12
            read (32) r4inp
            sstom12(:,:,it) = sst12(:,:,it)+r4inp
        end do

        if (iitest>=1) print *,' checking ocean model SST'
        call forchk (bmask_s,sstom12,ngp,12,100.0_dp,400.0_dp,273.0_dp)
    end if
end subroutine inbcon

subroutine forchk (fmask,field,ngp,nf,fmin,fmax,fset)
    ! Aux. routine forchk: Check consistency of sfc fields with land-sea mask
    ! and set undefined values to a constant (to avoid over/underflow)

    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: ngp, nf
    real(dp), intent(in) :: fmask(ngp)
    real(dp), intent(inout) :: field(ngp,nf)
    real(dp), intent(in) :: fmin, fmax, fset

    integer :: jf, jgp, nfault

    do jf = 1,nf
        nfault=0

        do jgp = 1,ngp
            if (fmask(jgp)>0.0_dp) then
                if (field(jgp,jf)<fmin .or. field(jgp,jf)>fmax) then
                    nfault = nfault+1
                end if
            else
                field(jgp,jf) = fset
            end if
        end do

        print *, ' field: ', jf, '   no. of faulty points:', nfault
    end do

    print *, ' undefined values set to', fset
end subroutine forchk

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
            if (itwn>itr) fsp(m,n)=zero
        end do
    end do

    call grid (fsp,fg2,1)
end subroutine truncg

subroutine fillsf(sf,ix,il,fmis)
    ! subroutine fillsf (sf,ix,nlat)
    ! Purpose: replace missing values in surface fields
    ! NB: it is assumed that non-missing values exist near the Equator

    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: ix, il
    real(dp), intent(inout) :: sf(ix,il)
    real(dp), intent(in) :: fmis

    real(dp) :: sf2(0:ix+1)
    integer :: khem, j, j1, j2, j3, i, nmis
    real(dp) :: fmean

    do khem = 1,2
       if (khem==1) then
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
                if (sf(i,j)<fmis) then
                    nmis = nmis+1
                    sf2(i) = 0.0_dp
                end if
            end do

            if (nmis<ix) fmean = sum(sf2(1:ix))/float(ix-nmis)

            do i=1,ix
                if (sf(i,j)<fmis) sf2(i) = fmean
            end do

            sf2(0)      = sf2(ix)
            sf2(ix+1) = sf2(1)
            do i=1,ix
                if (sf(i,j)<fmis) sf(i,j) = 0.5_dp*(sf2(i-1)+sf2(i+1))
            end do
        end do
    end do
end subroutine fillsf

subroutine NC_extract_variable(NCID, variable, it, inp)
    ! read a boundary field on from the opened netCDF file
    ! The input files are produces from GrADS files giving dimensions x,y,z,t
    ! The x and y dimensions match that of the run resolution, the z dimension
    ! is always 1 and the t dimension varies depending on the input.

    use netcdf
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: NCID, it
    character (len=*), intent(in) :: variable
    real(dp), intent(out) :: inp(ix,il,1,it)
    integer :: VARID

    call check(NF90_INQ_VARID(NCID, variable, VARID))
    call check(NF90_GET_VAR(NCID, VARID, inp))

    ! Fix undefined values
    where (inp <= -999) inp = 0.0_dp
end subroutine NC_extract_variable

subroutine check(status)
    ! Wrapper subroutine for NetCDF library functions that checks the returned
    ! status and stops the program on any errors
    use netcdf

    implicit none

    ! Status identifier output from netCDF library functions
    integer, intent (in) :: status

    if(status /= NF90_NOERR) then
        print *, NF90_STRERROR(status)
        stop "Stopped"
    end if
end subroutine check