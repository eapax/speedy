module mod_randfor
    use mod_atparam

    implicit none

    namelist /randfor/ nstrdf, indrdf

    ! Duration of random diabatic forcing ( 0 : no forcing, > 0 : no. of
    ! initial steps, < 0 : whole integration)
    integer :: nstrdf = 0

    ! Initialization index for random diabatic forcing
    integer :: indrdf = -1

    ! Flag for random diabatic forcing (set each timestep in dyn_stloop)
    logical :: lrandf

    ! Random diabatic forcing (initial. in INIRDF, modified by XS_RDF))
    real, dimension(:,:,:), allocatable :: randfh, randfv

    contains
        subroutine setup_randfor(fid)
            integer, intent(in) :: fid

            read(fid, randfor)
            write(*, randfor)

            allocate(randfh(ix,il,2))
            allocate(randfv(il,kx,2))
        end subroutine setup_randfor
end module
