!> @brief
!> Constants for large-scale condensation.
module mod_lsccon
    use mod_prec

    implicit none

    private
    public trlsc, rhlsc, drhlsc, rhblsc

    ! Relaxation time (in hours) for specific humidity 
    real(dp), parameter :: trlsc  = 4.0_dp

    ! Maximum relative humidity threshold (at sigma=1)
    real(dp), parameter :: rhlsc  = 0.9_dp

    ! Vertical range of relative humidity threshold
    real(dp), parameter :: drhlsc = 0.1_dp

    ! Relative humidity threshold for boundary layer
    real(dp), parameter :: rhblsc = 0.95_dp
end module
