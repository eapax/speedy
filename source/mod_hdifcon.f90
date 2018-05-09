module mod_hdifcon
    use mod_atparam
    use rp_emulator

    implicit none

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    type(rpe_var), allocatable, dimension(:,:) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    type(rpe_var), allocatable, dimension(:,:) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    type(rpe_var), allocatable, dimension(:) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    type(rpe_complex_var), allocatable, dimension(:,:) :: tcorh, qcorh

    contains
        subroutine setup_hdifcon()
            allocate(dmp(mx,nx))
            allocate(dmpd(mx,nx))
            allocate(dmps(mx,nx))
            allocate(dmp1(mx,nx))
            allocate(dmp1d(mx,nx))
            allocate(dmp1s(mx,nx))
            allocate(tcorv(kx))
            allocate(qcorv(kx))
            allocate(tcorh(mx,nx))
            allocate(qcorh(mx,nx))
        end subroutine setup_hdifcon

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
