subroutine ini_sea()
    use mod_cpl_flags, only: icsea
    use mod_atparam
    use mod_cli_sea, only: deglat_s
    use mod_var_sea
    use mod_prec, only: dp

    implicit none

    ! 1. Compute climatological fields for initial date
    call atm2sea(0)

    ! 2. Initialize prognostic variables of ocean/ice model
    !    in case of no restart or no coupling
    sst_om(:)  = sstcl_ob(:)      ! SST
    tice_om(:) = ticecl_ob(:)     ! sea ice temperature
    sice_om(:) = sicecl_ob(:)     ! sea ice fraction

    if (icsea<=0) sst_om(:) = 0.0_dp

    ! 3. Compute additional sea/ice variables
    wsst_ob(:) = 0.0_dp
    if (icsea>=4) call sea_domain('elnino',deglat_s,wsst_ob)

    call sea2atm(0)
end subroutine ini_sea

subroutine atm2sea(jday)
    ! subroutine atm2sea(jday)

    use mod_cpl_flags, only: icsea, icice, isstan
    use mod_atparam
    use mod_cplvar_sea, only: vsea_input
    use mod_date, only: iday, imont1, tmonth
    use mod_fluxes, only: hflux_s, hflux_i
    use mod_cli_sea, only: fmask_s, sst12, sice12, sstan3, hfseacl, sstom12
    use mod_var_sea, only: sstcl_ob, sicecl_ob, ticecl_ob, sstan_ob, sstcl_om,&
        & sst_om, tice_om
    use rp_emulator
    use mod_prec, only: dp,set_precision
 
    implicit none

    integer, intent(in) :: jday

    type(rpe_var) :: fmasks(ngp)        ! sea fraction
    type(rpe_var) :: hfyearm(ngp)       ! annual mean heat flux into the ocean
    integer :: j
    type(rpe_var) :: sstcl0, sstfr

    ! 1. Interpolate climatological fields and obs. SST anomaly
    !    to actual date

    print *, '--- Inside atm2sea'

    ! 1. Climatological SST
    call set_precision('rp_atm2sea_1') !And return it to 'normal'
    call forin5(ngp,imont1,tmonth,sst12,sstcl_ob)
    call set_precision('Default') !And return it to 'normal'


    ! 2. Climatological sea ice fraction
    call set_precision('rp_atm2sea_2') !And return it to 'normal'
    call forint(ngp,imont1,tmonth,sice12,sicecl_ob)
    call set_precision('Default') !And return it to 'normal'


    ! 3. SST anomaly
    call set_precision('rp_atm2sea_3') !And return it to 'normal'
    if (isstan>0) then
        if (iday==1 .and. jday>0) call OBS_SSTA
        call forint (ngp,2,tmonth,sstan3,sstan_ob)
    end if
    call set_precision('Default') !And return it to 'normal'


    ! 4. Ocean model climatological SST
    call set_precision('rp_atm2sea_4') !And return it to 'normal'
    if (icsea>=3) then
        call forin5 (ngp,imont1,tmonth,sstom12,sstcl_om)
    end if
    call set_precision('Default') !And return it to 'normal'


    ! 5. Adjust climatological fields over sea ice

    call set_precision('rp_atm2sea_5') !And return it to 'normal'
    ! SST at freezing point
    sstfr = rpe_literal(273.2_dp)-rpe_literal(1.8_dp)

    do j=1,ngp
        sstcl0 = sstcl_ob(j)

        if (sstcl_ob(j)>sstfr) then
            sicecl_ob(j) = min(rpe_literal(0.5_dp),sicecl_ob(j))
            ticecl_ob(j) = sstfr
            if (sicecl_ob(j)>rpe_literal(0.0_dp)) then
                sstcl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/&
                        (rpe_literal(1.0_dp)-sicecl_ob(j))
            end if
        else
            sicecl_ob(j) = max(rpe_literal(0.5_dp),sicecl_ob(j))
            ticecl_ob(j) = sstfr+(sstcl_ob(j)-sstfr)/sicecl_ob(j)
            sstcl_ob(j)  = sstfr
        end if

        if (icsea>=3) sstcl_om(j) = sstcl_om(j)+(sstcl_ob(j)-sstcl0)
    end do
    call set_precision('Default') !And return it to 'normal'

    !6. Reshaping
    call set_precision('rp_atm2sea_6') !And return it to 'normal'
    hfyearm = reshape(hfseacl, (/ngp/))
    fmasks = reshape(fmask_s, (/ngp/))
    call set_precision('Default') !And return it to 'normal'

    if (jday<=0) return
    
    call set_precision('rp_atm2sea_7') !And return it to 'normal'

        ! 2. Set input variables for mixed-layer/ocean model
        if (icsea>0 .or. icice>0) then
            vsea_input(:,1) = sst_om(:)
            vsea_input(:,2) = tice_om(:)
            vsea_input(:,3) = sicecl_ob(:)
            vsea_input(:,4) = hflux_s(:)
            vsea_input(:,5) = hflux_i(:)
            vsea_input(:,6) = sstcl_ob(:)
            vsea_input(:,7) = ticecl_ob(:)
            vsea_input(:,8) = hfyearm(:)
        end if
    call set_precision('Default') !And return it to 'normal'


    print *, 'Final outputs, summed, are:'
    print *, SUM(vsea_input(:,1)) 
    print *, SUM(vsea_input(:,2)) 
    print *, SUM(vsea_input(:,3)) 
    print *, SUM(vsea_input(:,4)) 
    print *, SUM(vsea_input(:,5)) 
    print *, SUM(vsea_input(:,6)) 
    print *, SUM(vsea_input(:,7))
    print *, SUM(vsea_input(:,8)) 
    print *, '------------'


        ! 3. Call message-passing routines to send data (if needed)
end subroutine atm2sea

subroutine sea2atm(jday)
    ! subroutine sea2atm(jday)

    use mod_cpl_flags, only: icsea, icice, isstan
    use mod_atparam
    use mod_cplvar_sea, only: vsea_output
    use mod_var_sea
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: jday

    if (jday>0 .and. (icsea>0 .or. icice>0)) then
        ! 1. Run ocean mixed layer or
        !    call message-passing routines to receive data from ocean model
        call sea_model

        ! 2. Get updated variables for mixed-layer/ocean model
        sst_om(:)   = vsea_output(:,1)      ! sst
        tice_om(:)  = vsea_output(:,2)      ! sea ice temperature
        sice_om(:)  = vsea_output(:,3)      ! sea ice fraction
    end if

    ! 3. Compute sea-sfc. anomalies and full fields for atm. model
    ! 3.1 SST
    sstan_am(:) = 0.0_dp

    if (icsea<=1) then
        if (isstan>0) sstan_am(:) = sstan_ob(:)

        ! Use observed SST (climatological or full field)
        sst_am(:) = sstcl_ob(:) + sstan_am(:)
    else if (icsea==2) then
        ! Use full ocean model SST
        sst_am(:) = sst_om(:)
    else if (icsea >= 3) then
        ! Define SST anomaly from ocean model ouput and climatology
        sstan_am(:) = sst_om(:) - sstcl_om(:)

        ! Merge with observed SST anomaly in selected area
        if (icsea>=4) then
            sstan_am(:) = sstan_am(:) + wsst_ob(:)*(sstan_ob(:)-sstan_am(:))
        end if

        ! Add observed SST climatology to model SST anomaly
        sst_am(:) = sstcl_ob(:) + sstan_am(:)
    end if

    ! 3.2 Sea ice fraction and temperature
    if (icice>0) then
        sice_am(:) = sice_om(:)
        tice_am(:) = tice_om(:)
    else
        sice_am(:) = sicecl_ob(:)
        tice_am(:) = ticecl_ob(:)
    end if

    sst_am(:)  = sst_am(:)+sice_am(:)*(tice_am(:)-sst_am(:))
    ssti_om(:) = sst_om(:)+sice_am(:)*(tice_am(:)-sst_om(:))
end subroutine sea2atm

subroutine rest_sea(imode)
    ! subroutine rest_sea(imode)

    ! Purpose : read/write sea variables from/to a restart file
    ! Input :   IMODE = 0 : read model variables from a restart file
    !                 = 1 : write model variables  to a restart file

    use mod_cpl_flags, only: icsea, icice
    use mod_atparam
    use mod_var_sea, only: sst_om, tice_om, sice_om, sst_am, tice_am, sice_am
    use mod_downscaling, only: ix_in, il_in, regrid
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    integer, intent(in) :: imode

    type(rpe_var) :: sst_c(ngp)              ! sst corrected for sea-ice values
    type(rpe_var) :: sstfr

    ! Sea variables at input resolution
    ! Data loaded in at full precision
    real(dp) :: sst_om_in(ix_in*il_in)
    real(dp) :: tice_om_in(ix_in*il_in)
    real(dp) :: sice_om_in(ix_in*il_in)

    if (imode==0) then
        ! Load data at full precision
        read (3)  sst_om_in     ! sst
        read (3) tice_om_in     ! sea ice temperature
        read (3) sice_om_in     ! sea ice fraction

        ! Interpolate to new grid
        if (ix_in/=ix .or. il_in/=il) then
            call regrid(sst_om_in, sst_om%val)
            call apply_truncation(sst_om)
        else
            sst_om = sst_om_in
        end if

        if (ix_in/=ix .or. il_in/=il) then
            call regrid(tice_om_in, tice_om%val)
            call apply_truncation(tice_om)
        else
            tice_om = tice_om_in
        end if

        if (ix_in/=ix .or. il_in/=il) then
            call regrid(sice_om_in, sice_om%val)
            call apply_truncation(sice_om)
        else
            sice_om = sice_om
        end if
    else
        !    write sea/ice model variables from coupled runs,
        !    otherwise write fields used by atmospheric model
        sstfr = 273.2_dp-1.8_dp

        if (icsea>0) then
            write (10) sst_om(:)%val
        else
            sst_c(:) = max(sst_am(:),sstfr)
            write (10) sst_c(:)%val
        end if

        if (icice>0) then
            write (10) tice_om(:)%val
            write (10) sice_om(:)%val
        else
            write (10) tice_am(:)%val
            write (10) sice_am(:)%val
        end if
    end if
end subroutine rest_sea

subroutine obs_ssta()
    ! subroutine obs_ssta

    ! Purpose : update observed SST anomaly array

    use netcdf
    use mod_atparam
    use mod_cli_sea, only: sstan3, bmask_s
    use mod_date, only: imonth, iyear, issty0
    use ppo_output_stream, only: check
    use mod_prec, only: dp
    use rp_emulator

    implicit none

    integer :: next_month
    integer :: NCID, VARID
    real(dp)   :: inp(ix,il)

    sstan3(:,:,1) = sstan3(:,:,2)
    sstan3(:,:,2) = sstan3(:,:,3)

    ! Compute next month given initial SST year
    next_month = (iyear - issty0) * 12 + imonth

    ! Read next month SST anomalies
    call check(NF90_OPEN('anomalies.nc', NF90_NOWRITE, NCID))
    call check(NF90_INQ_VARID(NCID, 'ssta', VARID))
    call check(NF90_GET_VAR(NCID, VARID, inp, &
                    start=(/1, 1, 1, 1+next_month/), count=(/ix, il, 1, 1/) ))
    call check(NF90_CLOSE(NCID))

    sstan3(1:ix,1:il,3)   = inp

    call forchk(bmask_s,sstan3(:,:,3),ngp,1,-50.0_dp,50.0_dp,0.0_dp)
end subroutine obs_ssta
