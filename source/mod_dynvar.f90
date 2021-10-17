!> @brief
!> Prognostic spectral variables for model dynamics, and geopotential.
!> Initialised in invars.
module mod_dynvar
    use mod_atparam
    use rp_emulator

    implicit none

    ! Prognostic spectral variables (updated in step)
    ! Vorticity
    type(rpe_complex_var), allocatable :: vor(:,:,:,:)

    ! Divergence
    type(rpe_complex_var), allocatable :: div(:,:,:,:)

    ! Absolute temperature
    type(rpe_complex_var), allocatable :: t(:,:,:,:)

    ! Copy of Absolute temperature used for debugging
    type(rpe_complex_var), allocatable :: tcopy(:,:,:,:)
    

    ! Log of (norm.) sfc pressure (p_s/p0)
    type(rpe_complex_var), allocatable :: PS(:,:,:)

    ! Tracers (tr.1: spec. humidity in g/kg)
    type(rpe_complex_var), allocatable :: TR(:,:,:,:,:)

    ! Geopotential (updated in geop)
    ! Atmos. geopotential
    type(rpe_complex_var), allocatable :: PHI(:,:,:)

    ! Surface geopotential
    type(rpe_complex_var), allocatable :: PHIS(:,:)



    contains
        subroutine setup_dynvar()
            allocate(vor(MX,NX,KX,2))
            allocate(div(MX,NX,KX,2))
            allocate(t(MX,NX,KX,2))
            allocate(tcopy(MX,NX,KX,2))
            allocate(PS(MX,NX,2))
            allocate(TR(MX,NX,KX,2,NTR))
            allocate(PHI(MX,NX,KX))
            allocate(PHIS(MX,NX))

        end subroutine setup_dynvar
end module
