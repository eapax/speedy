module mod_cplcon_sea
    use mod_atparam
    use rp_emulator
    use mod_prec

    implicit none

    namelist /sea/ depth_ml, dept0_ml, depth_ice, dept0_ice, &
            tdsst, tdice, fseamin, beta, &
            l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

    ! Constant parameters and fields in sea/ice model
    ! 1./heat_capacity (sea)
    type(rpe_var), allocatable :: rhcaps(:,:)

    ! 1./heat_capacity (ice)
    type(rpe_var), allocatable :: rhcapi(:,:)

    ! 1./dissip_time (sea)
    type(rpe_var), allocatable :: cdsea(:,:)

    ! 1./dissip_time (ice)
    type(rpe_var), allocatable :: cdice(:,:)

    ! ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
    type(rpe_var) :: depth_ml  ! High-latitude depth
    type(rpe_var) :: dept0_ml  ! Minimum depth (tropics)

    ! sea-ice depth : d + (d0-d)*(cos_lat)^2
    type(rpe_var) :: depth_ice ! High-latitude depth
    type(rpe_var) :: dept0_ice ! Minimum depth

    ! Dissipation time (days) for sea-surface temp. anomalies
    type(rpe_var) :: tdsst

    ! Dissipation time (days) for sea-ice temp. anomalies
    type(rpe_var) :: tdice

    ! Minimum fraction of sea for the definition of anomalies
    type(rpe_var) :: fseamin

    ! Heat flux coef. at sea/ice int.
    real(dp) :: beta = 1.0

    ! Geographical domain
    ! note : more than one regional domain may be set .true.
    logical :: l_globe  = .true.  ! global domain
    logical :: l_northe = .false. ! Northern hem. oceans (lat > 20N)
    logical :: l_natlan = .false. ! N. Atlantic (lat 20-80N, lon 100W-45E)
    logical :: l_npacif = .false. ! N. Pacific  (lat 20-80N, lon 100E-100W)
    logical :: l_tropic = .false. ! Tropics (lat 30S-30N)
    logical :: l_indian = .false. ! Indian Ocean (lat 30S-30N, lon 30-120E)

    contains
        subroutine setup_sea(fid)
            integer, intent(in) :: fid
            ! 1./heat_capacity (sea)
            allocate(rhcaps(ix,il))

            ! 1./heat_capacity (ice)
            allocate(rhcapi(ix,il))

            ! 1./dissip_time (sea)
            allocate(cdsea(ix,il))

            ! 1./dissip_time (ice)
            allocate(cdice(ix,il))

            read(fid, sea)
            fseamin = 1.0/fseamin

            write(*, sea)
        end subroutine setup_sea

        subroutine truncate_cplcon_sea()
            rhcaps = rhcaps
            rhcapi = rhcapi
            cdsea = cdsea
            cdice = cdice
        end subroutine truncate_cplcon_sea
end module
