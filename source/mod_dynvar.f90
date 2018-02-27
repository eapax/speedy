!> @brief
!> Prognostic spectral variables for model dynamics, and geopotential.
!> Initialised in invars.
module mod_dynvar
    use mod_atparam

    implicit none

    ! Prognostic spectral variables (updated in step)
    ! Vorticity
    complex, allocatable :: vor(:,:,:,:)

    ! Divergence 
    complex, allocatable :: div(:,:,:,:)

    ! Absolute temperature
    complex, allocatable :: t(:,:,:,:)

    ! Log of (norm.) sfc pressure (p_s/p0)
    complex, allocatable :: PS(:,:,:)

    ! Tracers (tr.1: spec. humidity in g/kg)
    complex, allocatable :: TR(:,:,:,:,:)

    ! Geopotential (updated in geop)
    ! Atmos. geopotential
    complex, allocatable :: PHI(:,:,:)

    ! Surface geopotential
    complex, allocatable :: PHIS(:,:)

    contains
        subroutine setup_dynvar()
            allocate(vor(MX,NX,KX,2))
            allocate(div(MX,NX,KX,2))
            allocate(t(MX,NX,KX,2))
            allocate(PS(MX,NX,2))
            allocate(TR(MX,NX,KX,2,NTR))
            allocate(PHI(MX,NX,KX))
            allocate(PHIS(MX,NX))
        end subroutine setup_dynvar
end module
