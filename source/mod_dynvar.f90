!> @brief
!> Prognostic spectral variables for model dynamics, and geopotential.
!> Initialised in invars.
module mod_dynvar
    use mod_atparam
    use mod_prec, only: dp

    implicit none

    ! Prognostic spectral variables (updated in step)
    ! Vorticity
    complex(dp), allocatable :: vor(:,:,:,:)

    ! Divergence
    complex(dp), allocatable :: div(:,:,:,:)

    ! Absolute temperature
    complex(dp), allocatable :: t(:,:,:,:)

    ! Log of (norm.) sfc pressure (p_s/p0)
    complex(dp), allocatable :: PS(:,:,:)

    ! Tracers (tr.1: spec. humidity in g/kg)
    complex(dp), allocatable :: TR(:,:,:,:,:)

    ! Geopotential (updated in geop)
    ! Atmos. geopotential
    complex(dp), allocatable :: PHI(:,:,:)

    ! Surface geopotential
    complex(dp), allocatable :: PHIS(:,:)

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
