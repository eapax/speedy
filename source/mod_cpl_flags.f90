!> @brief
!> Flags to set coupling options (see doc_instep.txt).
module mod_cpl_flags
    implicit none

    namelist /coupling_flags/ icland, icsea, icice, isstan

    ! Flag for land-coupling
    ! 0=no, 1=land-model
    integer :: icland = 1

    ! Flag for sea (SST) coupling
    ! 0 = precribed SST, no coupling
    ! 1 = precribed SST, ocean model forced by atm.
    ! 2 = full (uncorrected) SST from coupled ocean model
    ! 3 = SST anomaly from coupled ocean model + obs. SST clim.
    ! 4 = as 3 with prescribed SST anomaly in ElNino region
    integer :: icsea  = 0

    ! Flag for sea-ice coupling
    ! 0=no, 1=ice-model
    integer :: icice  = 1

    ! Flag for observed SST anomaly
    ! 0 = no (clim. SST), 1 = observed anomaly
    ! (active if ICSEA = 0, 1; set to 1 if ICSEA = 4)
    integer :: isstan = 1

    contains
        subroutine setup_coupling_flags(fid)
            integer, intent(in) :: fid

            read(fid, coupling_flags)

            write(*, coupling_flags)
        end subroutine setup_coupling_flags
end module
