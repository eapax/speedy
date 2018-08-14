module phy_cloud
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_cloud_parameters, ini_cloud, truncate_cloud, cloud
    public albcl, albcls, qcloud

    namelist /cloud_parameters/ &
            rhcl1, rhcl2, qacl, wpcl, pmaxcl, &
            clsmax, clsminl, gse_s0, gse_s1, &
            albcl, albcls

    ! rhcl1  = relative hum. threshold corr. to cloud cover = 0
    type(rpe_var) :: rhcl1
    ! rhcl2  = relative hum. corr. to cloud cover = 1
    type(rpe_var) :: rhcl2
    ! qacl   = specific hum. threshold for cloud cover
    type(rpe_var) :: qacl
    ! wpcl   = cloud c. weight for the sq. root of precip. (for p = 1 mm/day)
    type(rpe_var) :: wpcl
    ! pmaxcl = max. value of precip. (mm/day) contributing to cloud cover
    type(rpe_var) :: pmaxcl

    ! clsmax = maximum stratiform cloud cover
    type(rpe_var) :: clsmax
    ! clsminl= minimum stratiform cloud cover over land (for RH = 1)
    type(rpe_var) :: clsminl
    ! gse_s0 = gradient of dry static energy corresp. to strat.c.c. = 0
    type(rpe_var) :: gse_s0
    ! gse_s1 = gradient of dry static energy corresp. to strat.c.c. = 1
    type(rpe_var) :: gse_s1

    ! albcl  = cloud albedo (for cloud cover = 1)
    type(rpe_var) :: albcl
    ! albcls = stratiform cloud albedo (for st. cloud cover = 1)
    type(rpe_var) :: albcls

    ! Radiative properties of clouds (updated in cloud)
    ! qcloud = Equivalent specific humidity of clouds
    type(rpe_var), dimension(:), allocatable :: qcloud

    ! Local derived variables
    type(rpe_var) :: rrcl, albcor, clfact

    contains
        subroutine setup_cloud_parameters(fid)
            ! Read namelist variables
            integer, intent(in) :: fid

            read(fid, cloud_parameters)
            write(*, cloud_parameters)

            allocate(qcloud(ngp))
        end subroutine setup_cloud_parameters

        subroutine ini_cloud()
            rrcl = 1.0_dp/(rhcl2-rhcl1)
            albcor  = albcl/0.5_dp
            clfact = 1.2_dp

        end subroutine ini_cloud

        subroutine truncate_cloud()
            ! Truncate local variables for cloud scheme
            ! Namelist variables
            call apply_truncation(rhcl1)
            call apply_truncation(rhcl2)
            call apply_truncation(qacl)
            call apply_truncation(wpcl)
            call apply_truncation(pmaxcl)
            call apply_truncation(clsmax)
            call apply_truncation(clsminl)
            call apply_truncation(gse_s0)
            call apply_truncation(gse_s1)
            call apply_truncation(albcl)
            call apply_truncation(albcls)

            ! Local derived variables
            call apply_truncation(rrcl)
            call apply_truncation(albcor)
            call apply_truncation(clfact)
        end subroutine truncate_cloud

        subroutine cloud(qa,rh,precnv,precls,iptop,gse,fmask,icltop,cloudc,clstr)
            !  subroutine cloud (qa,rh,precnv,precls,iptop,gse,fmask,
            ! &                  icltop,cloudc,clstr)
            !
            !  Purpose: Compute cloud-top level and cloud cover
            !  Input:   qa     = specific humidity [g/kg]                (3-dim)
            !           rh     = relative humidity                       (3-dim)
            !           precnv = convective precipitation                (2-dim)
            !           precls = large-scale precipitation               (2-dim)
            !           iptop  = top level of precipitating cloud        (2-dim)
            !           gse    = gradient of dry st. energy (dSE/dPHI)   (2-dim)
            !           fmask  = fractional land-sea mask                (2-dim)
            !  Output:  icltop = cloud top level (all clouds)            (2-dim)
            !           cloudc = total cloud cover                       (2-dim)
            !           clstr  = stratiform cloud cover                  (2-dim)

            integer :: iptop(ngp)
            type(rpe_var), intent(in) :: qa(ngp,kx), rh(ngp,kx), precnv(ngp), &
                    precls(ngp), gse(ngp), fmask(ngp)
            type(rpe_var), intent(out) :: cloudc(ngp), clstr(ngp)
            integer, intent(out) :: icltop(ngp)

            integer :: inew, j, k
            type(rpe_var) :: cl1, clstrl, drh, fstab, pr1, rgse

            ! 1.  Cloud cover, defined as the sum of:
            !     - a term proportional to the square-root of precip. rate
            !     - a quadratic function of the max. relative humidity
            !       in tropospheric layers above PBL where Q > QACL :
            !       ( = 0 for RHmax < RHCL1, = 1 for RHmax > RHCL2 )
            !     Cloud-top level: defined as the highest (i.e. least sigma)
            !       between the top of convection/condensation and
            !       the level of maximum relative humidity.

            do j=1,ngp
                if (rh(j,kxm) > rhcl1) then
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
                    if (drh.gt.cloudc(j).and.qa(j,k).gt.qacl) then
                        cloudc(j) = drh
                        icltop(j) = k
                    end if
                end do
            end do

            do j=1,ngp
                cl1 = min(rpe_literal(1.0_dp),cloudc(j)*rrcl)
                pr1 = min(pmaxcl,rpe_literal(86.4_dp)*(precnv(j)+precls(j)))
                cloudc(j) = min(rpe_literal(1.0_dp),wpcl*sqrt(pr1)+cl1*cl1)
                icltop(j) = min(iptop(j),icltop(j))
            end do

            ! 2.  Equivalent specific humidity of clouds
            qcloud = qa(:,kxm)

            ! 3. Stratiform clouds at the top of PBL
            inew = 1

            if (inew.gt.0) then
                !        CLSMAX  = 0.6
                !        CLSMINL = 0.15
                !        GSE_S0  = 0.25
                !        GSE_S1  = 0.40


                rgse   = rpe_literal(1.0_dp)/(gse_s1-gse_s0)

                do j=1,ngp
                    ! Stratocumulus clouds over sea
                    fstab    = max(rpe_literal(0.0_dp), &
                            min(rpe_literal(1.0_dp), rgse*(gse(j)-gse_s0)))
                    clstr(j) = fstab*&
                            max(clsmax-clfact*cloudc(j), rpe_literal(0.0_dp))
                    ! Stratocumulus clouds over land
                    clstrl   = max(clstr(j),clsminl)*rh(j,kx)
                    clstr(j) = clstr(j)+fmask(j)*(clstrl-clstr(j))
                end do
            else
                clsmax  = 0.3_dp
                clsminl = 0.1_dp


                do j=1,ngp
                    ! stratocumulus clouds over sea
                    clstr(j) = max(clsmax-cloudc(j),rpe_literal(0.0_dp))
                    ! rescale for consistency with previous albedo values
                    clstr(j) = clstr(j)*albcor
                    ! correction for aerosols over land
                    clstr(j) = clstr(j)+fmask(j)*(clsminl-clstr(j))
                end do
            end if
        end subroutine cloud
end module phy_cloud
