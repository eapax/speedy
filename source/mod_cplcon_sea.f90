module mod_cplcon_sea
    use mod_atparam

    implicit none

    private
    public rhcaps, rhcapi, cdsea, cdice, &
            depth_ml, dept0_ml, depth_ice, dept0_ice, &
            tdsst, tdice, fseamin, beta, &
            l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

    namelist /sea_model/ depth_ml, dept0_ml, depth_ice, dept0_ice, &
            tdsst, tdice, fseamin, beta, &
            l_globe, l_northe, l_natlan, l_npacif, l_tropic, l_indian

    ! Constant parameters and fields in sea/ice model
    ! 1./heat_capacity (sea)
    real :: rhcaps(ix,il)

    ! 1./heat_capacity (ice)
    real :: rhcapi(ix,il)

    ! 1./dissip_time (sea)
    real :: cdsea(ix,il)

    ! 1./dissip_time (ice)
    real :: cdice(ix,il)

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
    logical, parameter :: l_globe  = .true.  ! global domain
    logical, parameter :: l_northe = .false. ! Northern hem. oceans (lat > 20N)
    logical, parameter :: l_natlan = .false. ! N. Atlantic (lat 20-80N, lon 100W-45E)
    logical, parameter :: l_npacif = .false. ! N. Pacific  (lat 20-80N, lon 100E-100W)
    logical, parameter :: l_tropic = .false. ! Tropics (lat 30S-30N)
    logical, parameter :: l_indian = .false. ! Indian Ocean (lat 30S-30N, lon 30-120E)

    contains
        subroutine setup_sea_model(fid)
            integer, intent(in) :: fid

            read(fid, sea_model)
            fseamin = 1.0/fseamin
        end subroutine setup_sea_model
end module
