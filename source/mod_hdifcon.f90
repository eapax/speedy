module mod_hdifcon
    use mod_atparam

    implicit none

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    real, allocatable, dimension(:,:) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    real, allocatable, dimension(:,:) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    real, allocatable, dimension(:) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    complex, allocatable, dimension(:,:) :: tcorh, qcorh

    contains
        subroutine setup_hdifcon()
            ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
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
end module
