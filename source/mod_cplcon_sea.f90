module mod_cplcon_sea
    use mod_atparam

    implicit none

    namelist /sea/ depth_ml, dept0_ml, depth_ice, dept0_ice, &
            tdsst, tdice, fseamin, beta, &
            l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

    ! Constant parameters and fields in sea/ice model
    ! 1./heat_capacity (sea)
    real, allocatable :: rhcaps(:,:)

    ! 1./heat_capacity (ice)
    real, allocatable :: rhcapi(:,:)

    ! 1./dissip_time (sea)
    real, allocatable :: cdsea(:,:)

    ! 1./dissip_time (ice)
    real, allocatable :: cdice(:,:)

    ! ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
    real :: depth_ml = 60.               ! High-latitude depth
    real :: dept0_ml = 40.               ! Minimum depth (tropics)

    ! sea-ice depth : d + (d0-d)*(cos_lat)^2
    real :: depth_ice = 2.5              ! High-latitude depth
    real :: dept0_ice = 1.5              ! Minimum depth

    ! Dissipation time (days) for sea-surface temp. anomalies
    real :: tdsst  = 90.

    ! Dissipation time (days) for sea-ice temp. anomalies
    real :: tdice = 30.

    ! Minimum fraction of sea for the definition of anomalies
    real :: fseamin = 3.0

    ! Heat flux coef. at sea/ice int.
    real :: beta = 1.0

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
end module
