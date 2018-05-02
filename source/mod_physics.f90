!> @brief
!> Logical flags to control certain behaviour.
module mod_physics
    implicit none

    namelist /physics/ lco2, sppt_on, nstrad
    ! Logical flags to activate processes throughout the integration
    ! Flag for CO2 optical thickness increase
    logical :: lco2 = .false.

    ! Turn on SPPT?
    logical :: sppt_on = .false.

    ! Period (no. of steps) for shortwave radiation
    integer :: nstrad = 3

    ! Logical flags to activate processes in selected time steps (updated in
    ! STLOOP)
    ! Flag for shortwave radiation routine
    logical :: lradsw

    contains
        subroutine setup_physics(fid)
            integer, intent(in) :: fid

            read(fid, physics)

            write(*, physics)
        end subroutine setup_physics
end module
