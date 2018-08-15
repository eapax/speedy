module phy_cloud
    use mod_atparam
    use rp_emulator
    use mod_prec, only: dp

    implicit none

    private
    public setup_cloud_parameters, ini_cloud, truncate_cloud, cloud

    namelist /cloud_parameters/ &
            rhcl1, rhcl2, qacl, wpcl, pmaxcl, &
            clsmax, clsminl, gse_s0, gse_s1

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

    ! Local derived variables
    type(rpe_var) :: rrcl, clfact

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

            ! Local derived variables
            call apply_truncation(rrcl)
            call apply_truncation(clfact)
        end subroutine truncate_cloud

        subroutine cloud(&
                qa_in, rh_in, precnv_in, precls_in, iptop, gse_in, fmask_in,&
                icltop, cloudc, clstr)
            !  subroutine cloud (qa,rh,precnv,precls,iptop,gse,fmask,
            ! &                  icltop,cloudc,clstr)
            !
            !  Purpose: Compute cloud-top level and cloud cover
            !  Input:   qa     = specific humidity [g/kg]                (3-dim)
            type(rpe_var), intent(in) :: qa_in(ngp,kx)
            !           rh     = relative humidity                       (3-dim)
            type(rpe_var), intent(in) :: rh_in(ngp,kx)
            !           precnv = convective precipitation                (2-dim)
            type(rpe_var), intent(in) :: precnv_in(ngp)
            !           precls = large-scale precipitation               (2-dim)
            type(rpe_var), intent(in) :: precls_in(ngp)
            !           iptop  = top level of precipitating cloud        (2-dim)
            integer, intent(in) :: iptop(ngp)
            !           gse    = gradient of dry st. energy (dSE/dPHI)   (2-dim)
            type(rpe_var), intent(in) :: gse_in(ngp)
            !           fmask  = fractional land-sea mask                (2-dim)
            type(rpe_var), intent(in) :: fmask_in(ngp)
            !  Output:  icltop = cloud top level (all clouds)            (2-dim)
            integer, intent(out) :: icltop(ngp)
            !           cloudc = total cloud cover                       (2-dim)
            type(rpe_var), intent(out) :: cloudc(ngp)
            !           clstr  = stratiform cloud cover                  (2-dim)
            type(rpe_var), intent(out) :: clstr(ngp)

            ! Local copies of input variables
            type(rpe_var) :: qa(ngp,kx), rh(ngp,kx), &
                    precnv(ngp), precls(ngp), gse(ngp), fmask(ngp)

            ! Local variables
            integer :: j, k
            type(rpe_var) :: cl1, clstrl, drh, fstab, pr1, rgse

            ! 0. Pass input variables to local copies, triggering call to
            !    apply_truncation
            qa = qa_in
            rh = rh_in
            precnv = precnv_in
            precls = precls_in
            gse = gse_in
            fmask = fmask_in

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
                    if (drh > cloudc(j).and.qa(j,k) > qacl) then
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

            ! 3. Stratiform clouds at the top of PBL
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
        end subroutine cloud
end module phy_cloud
