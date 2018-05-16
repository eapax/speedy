subroutine dmflux(iadd)
    ! subroutine dmflux (iadd)
    !
    ! Purpose: Add up fluxes to provide daily averages 
    !          used in sea/land models and daily/time-mean output
    ! Input: IADD = 0 to initialize storage arrays to 0
    !             > 0 to increment arrays with current flux values  

    use mod_tsteps, only: nsteps
    use mod_atparam
    use mod_flx_land
    use mod_flx_sea
    use mod_physcon, only: alhc, sbc
    use mod_surfcon, only: fmask1
    use mod_var_sea, only: tice_am, sice_am
    use mod_physvar
    use phy_radiat, only: albsea, albice, emisfc
    use mod_date, only: ihour

    implicit none

    integer, intent(in) :: iadd
    integer :: j

    real :: difice(ngp)

    real :: fland(ngp), esbc, rsteps, sstfr, sstfr4

    fland = reshape(fmask1,(/ngp/))

    ! 1. Initialization
    if (iadd.le.0) then
        if (ihour /= 0) then
            ! Read from flux file
            open (100,file='fluxes.grd',form='unformatted',access='direct',&
                    recl=8*ngp)
            read (100,rec=4) (hflux_l(j),j=1,ngp)
            read (100,rec=14) (hflux_s(j),j=1,ngp)
            read (100,rec=15) (hflux_i(j),j=1,ngp)
            close (100)
        else
            ! Set all daily-mean arrays to zero
            hflux_l(:) = 0.
            hflux_s(:) = 0.
            hflux_i(:) = 0.
        end if
        return
    end if

    rsteps = 1./real(nsteps)

    ! SST at freezing point
    sstfr  = 273.2-1.8

    sstfr4 = sstfr**4
    esbc   = emisfc*sbc

    ! 2. Store fluxes over land (SI units, all heat fluxes downw.)
    hflux_l(:) = hflux_l(:) + hfluxn(:,1)*rsteps

    ! 3. Store fluxes over sea (SI units, all heat fluxes downw.)
    ! Difference in net (downw.) heat flux between ice and sea surface
    difice(:) = (albsea-albice)*ssrd(:)+ esbc*(sstfr4-tice_am(:)**4)&
        & + shf(:,2)+evap(:,2)*alhc

    hflux_s(:) = hflux_s(:) + rsteps* hfluxn(:,2)
    hflux_i(:) = hflux_i(:) + rsteps*(hfluxn(:,2)+difice(:)*(1.-sice_am(:)))

    ! 4.1 Store fluxes for daily-mean output

    ! Multiply net heat fluxes by land or sea fractions
    hfluxn(:,1) = hfluxn(:,1)*fland(:)
    hfluxn(:,2) = hfluxn(:,2)*(1.-fland(:))

    ! End of flux increment
end
