subroutine sea_model_init(fmask_s,rlat)
    !  subroutine sea_model_init (fmask_s,rlat)
    !
    !  Purpose : Initialization of sea model
    !  Initialized common blocks: sea_mc

    use mod_atparam
    use mod_cplcon_sea
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    ! Input variables
    type(rpe_var) :: fmask_s(ix,il)            ! sea mask (fraction of sea)
    type(rpe_var) :: rlat(il)                  ! latitudes in degrees

    ! Auxiliary variables

    ! Domain mask
    real(dp) :: dmask(ix,il)

    ! Heat capacity of mixed-l
    real(dp) :: hcaps(il)

    ! Heat capacity of sea-ice
    real(dp) :: hcapi(il)

    integer :: i, j
    real(dp) :: coslat, crad

    ! 1. Set geographical domain, heat capacities and dissipation times
    !    for sea (mixed layer) and sea-ice

    ! Heat capacities per m^2 (depth*heat_cap/m^3)
    crad=asin(1.0_dp)/90.0_dp
    do j=1,il
        coslat   = cos(crad*rlat(j))
        hcaps(j) = 4.18d+6*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
        hcapi(j) = 1.93d+6*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
    end do

    ! 3. Compute constant parameters and fields

    ! Set domain mask
    if (l_globe) then
        dmask(:,:) = 1.0_dp
    else
        dmask(:,:) = 0.0_dp
        if (l_northe) call SEA_DOMAIN ('northe',rlat,dmask)
        !fkif (l_arctic) call SEA_DOMAIN ('arctic',rlat,dmask)
        if (l_natlan) call SEA_DOMAIN ('natlan',rlat,dmask)
        if (l_npacif) call SEA_DOMAIN ('npacif',rlat,dmask)
        if (l_tropic) call SEA_DOMAIN ('tropic',rlat,dmask)
        if (l_indian) call SEA_DOMAIN ('indian',rlat,dmask)
    end if

    ! Smooth latitudinal boundaries and blank out land points
    do j=2,il-1
        rhcaps(:,j) = 0.25_dp*(dmask(:,j-1)+2*dmask(:,j)+dmask(:,j+1))
    end do
    dmask(:,2:il-1) = rhcaps(:,2:il-1)

    do j=1,il
        do i=1,ix
            if (fmask_s(i,j).lt.fseamin) dmask(i,j) = 0
        end do
    end do

    ! Set heat capacity and dissipation time over selected domain
    do j=1,il
        rhcaps(:,j) = 86400.0_dp/hcaps(j)
        rhcapi(:,j) = 86400.0_dp/hcapi(j)
    end do

    cdsea = dmask*tdsst/(1.0_dp+dmask*tdsst)
    cdice = dmask*tdice/(1.0_dp+dmask*tdice)
end

subroutine sea_model()
    ! subroutine sea_model

    ! Purpose : Integrate slab ocean and sea-ice models for one day

    use mod_atparam
    use mod_cplcon_sea
    use mod_cplvar_sea
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    ! Input variables:
    type(rpe_var) ::  sst0(ix,il)     ! SST at initial time
    type(rpe_var) :: tice0(ix,il)     ! sea ice temp. at initial time
    type(rpe_var) :: sice0(ix,il)     ! sea ice fraction at initial time
    type(rpe_var) :: hfsea(ix,il)     ! sea+ice  sfc. heat flux between t0 and t1
    type(rpe_var) :: hfice(ix,il)     ! ice-only sfc. heat flux between t0 and t1

    type(rpe_var) ::  sstcl1(ix,il)   ! clim. SST at final time
    type(rpe_var) :: ticecl1(ix,il)   ! clim. sea ice temp. at final time
    type(rpe_var) :: hfseacl(ix,il)   ! clim. heat flux due to advection/upwelling

    ! Output variables
    type(rpe_var) ::  sst1(ix,il)     ! SST at final time
    type(rpe_var) :: tice1(ix,il)     ! sea ice temp. at final time
    type(rpe_var) :: sice1(ix,il)     ! sea ice fraction at final time

    ! Auxiliary variables
    type(rpe_var) :: hflux(ix,il)   ! net sfc. heat flux
    type(rpe_var) :: tanom(ix,il)   ! sfc. temperature anomaly
    type(rpe_var) :: cdis(ix,il)    ! dissipation ceofficient

    type(rpe_var) :: anom0, sstfr

    sst0 = reshape(vsea_input(:,1), (/ix, il/))
    tice0 = reshape(vsea_input(:,2), (/ix, il/))
    sice0 = reshape(vsea_input(:,3), (/ix, il/))
    hfsea = reshape(vsea_input(:,4), (/ix, il/))
    hfice = reshape(vsea_input(:,5), (/ix, il/))
    sstcl1 = reshape(vsea_input(:,6), (/ix, il/))
    ticecl1 = reshape(vsea_input(:,7), (/ix, il/))
    hfseacl = reshape(vsea_input(:,8), (/ix, il/))

    sstfr = rpe_literal(273.2_dp)-rpe_literal(1.8_dp)    ! SST at freezing point

    ! 1. Ocean mixed layer
    ! Net heat flux
    hflux = hfsea-hfseacl-sice0*(hfice+beta*(sstfr-tice0))

    ! Anomaly at t0 minus climatological temp. tendency
    tanom = sst0 - sstcl1

    ! Time evoloution of temp. anomaly
    tanom = cdsea*(tanom+rhcaps*hflux)

    ! Full SST at final time
    sst1 = tanom + sstcl1

    ! 2. Sea-ice slab model

    ! Net heat flux
    hflux = hfice + beta*(sstfr-tice0)

    ! Anomaly w.r.t final-time climatological temp.
    tanom = tice0 - ticecl1

    ! Definition of non-linear damping coefficient
    anom0     = 20.0_dp
    cdis = cdice*(anom0/(anom0+abs(tanom)))
    !cdis(:,:) = cdice(:,:)

    ! Time evoloution of temp. anomaly
    tanom = cdis*(tanom+rhcapi*hflux)

    ! Full ice temperature at final time
    tice1 = tanom + ticecl1

    ! Persistence of sea ice fraction
    sice1 = sice0

    vsea_output(:,1) = reshape(sst1, (/ngp/))
    vsea_output(:,2) = reshape(tice1, (/ngp/))
    vsea_output(:,3) = reshape(sice1, (/ngp/))
end

subroutine sea_domain(cdomain,rlat,dmask)
    ! subroutine sea_domain (cdomain,rlat,dmask)

    ! Purpose : Definition of ocean domains

    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    ! Input variables

    character(len=6), intent(in) :: cdomain           ! domain name
    type(rpe_var), intent(in) :: rlat(il)               ! latitudes in degrees

    ! Output variables (initialized by calling routine)
    real(dp), intent(inout) :: dmask(ix,il)         ! domain mask

    integer :: i, j
    real(dp) :: arlat, dlon, rlon, rlonw, wlat

    print *, 'sea domain : ', cdomain

    dlon = rpe_literal(360.0_dp)/rpe_literal(ix)

    if (cdomain.eq.'northe') then
        do j=1,il
            if (rlat(j).gt.rpe_literal(20.0_dp)) dmask(:,j) = 1.0_dp
        end do
    end if

    if (cdomain.eq.'natlan') then
         do j=1,il
           if (rlat(j).gt.rpe_literal(20.0_dp) &
                   .and.rlat(j).lt.rpe_literal(80.0_dp)) then
             do i=1,ix
               rlon = (i-1)*dlon
               if (rlon.lt.rpe_literal(45.0_dp) &
                   .or.rlon.gt.rpe_literal(260.0_dp)) dmask(i,j) = 1.0_dp
             end do
           end if
         end do
    end if

    if (cdomain.eq.'npacif') then
        do j=1,il
            if (rlat(j).gt.rpe_literal(20.0_dp) &
                    .and.rlat(j).lt.rpe_literal(65.0_dp)) then
                do i=1,ix
                    rlon = (i-1)*dlon
                    if (rlon.gt.rpe_literal(120.0_dp) &
                        .and.rlon.lt.rpe_literal(260.0_dp)) dmask(i,j) = 1.0_dp
                end do
            end if
        end do
    end if

    if (cdomain.eq.'tropic') then
        do j=1,il
            if (rlat(j).gt.rpe_literal(-30.0_dp) &
                .and.rlat(j).lt.rpe_literal(30.0_dp)) dmask(:,j) = 1.0_dp
        end do
    end if

    if (cdomain.eq.'indian') then
        do j=1,il
            if (rlat(j).gt.rpe_literal(-30.0_dp) &
                    .and.rlat(j).lt.rpe_literal(30.0_dp)) then
                do i=1,ix
                    rlon = (i-1)*dlon
                    if (rlon.gt.rpe_literal(30.0_dp)&
                        .and.rlon.lt.rpe_literal(120.0_dp)) dmask(i,j) = 1.0_dp
                end do
            end if
        end do
    end if

    if (cdomain.eq.'elnino') then
        do j=1,il
            arlat = abs(rlat(j))
            if (arlat.lt.rpe_literal(25.0_dp)) then
                wlat = 1.0_dp
                if (arlat.gt.rpe_literal(15.0_dp)) then
                    wlat = (rpe_literal(0.1_dp)*(rpe_literal(25.0_dp)-arlat))**2
                end if
                rlonw = rpe_literal(300.0_dp)-2*max(rlat(j),rpe_literal(0.0_dp))
                do i=1,ix
                    rlon = (i-1)*dlon
                    if (rlon.gt.rpe_literal(165.0_dp).and.rlon.lt.rlonw) then
                        dmask(i,j) = wlat
                    else if (rlon.gt.rpe_literal(155.0_dp)&
                            .and.rlon.lt.rpe_literal(165.0_dp)) then
                        dmask(i,j) = wlat*rpe_literal(0.1_dp)*&
                                (rlon-rpe_literal(155.0_dp))
                    end if
                end do
            end if
        end do
    end if
end
