module mod_hdifcon
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Damping coef. for horizontal diffusion (explicit) (initial. in indyns)
    real(dp), allocatable, dimension(:,:) :: dmp, dmpd, dmps

    ! Damping coef. for horizontal diffusion (implicit) (initial. in indyns)
    real(dp), allocatable, dimension(:,:) :: dmp1, dmp1d, dmp1s

    ! Vertical comp. of orographic correction (initial. in INDYNS)
    real(dp), allocatable, dimension(:) :: tcorv, qcorv

    ! Horizontal component of orographic correction (updated in FORDATE)
    complex(dp), allocatable, dimension(:,:) :: tcorh, qcorh

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
end module
