module mod_randfor
    use mod_atparam

    implicit none

    ! Random diabatic forcing (initial. in INIRDF, modified by XS_RDF))
    real, dimension(:,:,:), allocatable :: randfh, randfv

    contains
        subroutine setup_randfor()
            allocate(randfh(ix,il,2))
            allocate(randfv(il,kx,2))
        end subroutine setup_randfor
end module
