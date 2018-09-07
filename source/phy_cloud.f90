module phy_cloud
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    private
    public setup_cloud_parameters, ini_cloud, cloud

    namelist /cloud_parameters/ &
            rhcl1, rhcl2, qacl, wpcl, pmaxcl, &
            clsmax, clsminl, gse_s0, gse_s1

    ! rhcl1  = relative hum. threshold corr. to cloud cover = 0
    real(dp) :: rhcl1
    ! rhcl2  = relative hum. corr. to cloud cover = 1
    real(dp) :: rhcl2
    ! qacl   = specific hum. threshold for cloud cover
    real(dp) :: qacl
    ! wpcl   = cloud c. weight for the sq. root of precip. (for p = 1 mm/day)
    real(dp) :: wpcl
    ! pmaxcl = max. value of precip. (mm/day) contributing to cloud cover
    real(dp) :: pmaxcl

    ! clsmax = maximum stratiform cloud cover
    real(dp) :: clsmax
    ! clsminl= minimum stratiform cloud cover over land (for RH = 1)
    real(dp) :: clsminl
    ! gse_s0 = gradient of dry static energy corresp. to strat.c.c. = 0
    real(dp) :: gse_s0
    ! gse_s1 = gradient of dry static energy corresp. to strat.c.c. = 1
    real(dp) :: gse_s1

    ! Local derived variables
    real(dp) :: rrcl, clfact

    contains
        subroutine setup_cloud_parameters(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, cloud_parameters)
            write(*, cloud_parameters)
        end subroutine setup_cloud_parameters

        subroutine ini_cloud()
            rrcl = 1.0_dp/(rhcl2-rhcl1)
            clfact = 1.2_dp
        end subroutine ini_cloud

        subroutine cloud(&
                qa, rh, precnv, precls, iptop, gse, fmask,&
                icltop, cloudc, clstr)
            !  subroutine cloud (qa,rh,precnv,precls,iptop,gse,fmask,
            ! &                  icltop,cloudc,clstr)
            !
            !  Purpose: Compute cloud-top level and cloud cover
            !  Input:   qa     = specific humidity [g/kg]                (3-dim)
            real(dp), intent(in) :: qa(ngp,kx)
            !           rh     = relative humidity                       (3-dim)
            real(dp), intent(in) :: rh(ngp,kx)
            !           precnv = convective precipitation                (2-dim)
            real(dp), intent(in) :: precnv(ngp)
            !           precls = large-scale precipitation               (2-dim)
            real(dp), intent(in) :: precls(ngp)
            !           iptop  = top level of precipitating cloud        (2-dim)
            integer, intent(in) :: iptop(ngp)
            !           gse    = gradient of dry st. energy (dSE/dPHI)   (2-dim)
            real(dp), intent(in) :: gse(ngp)
            !           fmask  = fractional land-sea mask                (2-dim)
            real(dp), intent(in) :: fmask(ngp)
            !  Output:  icltop = cloud top level (all clouds)            (2-dim)
            integer, intent(out) :: icltop(ngp)
            !           cloudc = total cloud cover                       (2-dim)
            real(dp), intent(out) :: cloudc(ngp)
            !           clstr  = stratiform cloud cover                  (2-dim)
            real(dp), intent(out) :: clstr(ngp)

            ! Local variables
            integer :: j, k
            real(dp) :: cl1, clstrl, drh, fstab, pr1, rgse

            ! 1.  Cloud cover, defined as the sum of:
            !     - a term proportional to the square-root of precip. rate
            !     - a quadratic function of the max. relative humidity
            !       in tropospheric layers above PBL where Q > QACL :
            !       ( = 0 for RHmax < RHCL1, = 1 for RHmax > RHCL2 )
            !     Cloud-top level: defined as the highest (i.e. least sigma)
            !       between the top of convection/condensation and
            !       the level of maximum relative humidity.

            do j=1,ngp
                if (rh(j,kxm)>rhcl1) then
                    cloudc(j) = rh(j,kxm)-rhcl1
                    icltop(j) = kxm
                else
                    cloudc(j) = 0.0_dp
                    icltop(j) = kxp
                end if
            end do

            do k=3,kx-2
                do j=1,ngp
                    drh = rh(j,k)-rhcl1
                    if (drh>cloudc(j) .and. qa(j,k)>qacl) then
                        cloudc(j) = drh
                        icltop(j) = k
                    end if
                end do
            end do

            do j=1,ngp
                cl1 = min(1.0_dp,cloudc(j)*rrcl)
                pr1 = min(pmaxcl,86.4_dp*(precnv(j)+precls(j)))
                cloudc(j) = min(1.0_dp,wpcl*sqrt(pr1)+cl1*cl1)
                icltop(j) = min(iptop(j),icltop(j))
            end do

            ! 3. Stratiform clouds at the top of PBL
            rgse   = 1.0_dp/(gse_s1-gse_s0)

            do j=1,ngp
                ! Stratocumulus clouds over sea
                fstab    = max(0.0_dp, &
                        min(1.0_dp, rgse*(gse(j)-gse_s0)))
                clstr(j) = fstab*&
                        max(clsmax-clfact*cloudc(j), 0.0_dp)
                ! Stratocumulus clouds over land
                clstrl   = max(clstr(j),clsminl)*rh(j,kx)
                clstr(j) = clstr(j)+fmask(j)*(clstrl-clstr(j))
            end do
        end subroutine cloud
end module phy_cloud
