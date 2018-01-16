module mod_hdifcon
    use mod_atparam
    use rp_emulator

    implicit none

    private
    public truncate_hdifcon
    public dmp, dmpd, dmps, dmp1, dmp1d, dmp1s, tcorv, qcorv, tcorh, qcorh

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    type(rpe_var), dimension(mx,nx) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    type(rpe_var), dimension(mx,nx) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    type(rpe_var), dimension(kx) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    type(rpe_complex_var), dimension(mx,nx) :: tcorh, qcorh

    contains
        subroutine truncate_hdifcon()
            dmp = dmp
            dmpd = dmpd
            dmps = dmps
            dmp1 = dmp1
            dmp1d = dmp1d
            dmp1s = dmp1s
            tcorv = tcorv
            qcorv = qcorv
            tcorh = tcorh
            qcorh = qcorh
        end subroutine truncate_hdifcon
end module
