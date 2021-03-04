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

    ! Namelist parameters used to set up model constants
    ! ocean mixed layer depth: d + (d0-d)*(cos_lat)^3
    real(dp) :: depth_ml  ! High-latitude depth
    real(dp) :: dept0_ml  ! Minimum depth (tropics)
    ! sea-ice depth : d + (d0-d)*(cos_lat)^2
    real(dp) :: depth_ice ! High-latitude depth
    real(dp) :: dept0_ice ! Minimum depth
    ! Dissipation time (days) for sea-surface temp. anomalies
    real(dp) :: tdsst
    ! Dissipation time (days) for sea-ice temp. anomalies
    real(dp) :: tdice
    ! Minimum fraction of sea for the definition of anomalies
    real(dp) :: fseamin

    ! Namelist parameters used to integrate sea model
    ! Heat flux coef. at sea/ice int.
    type(rpe_var) :: beta

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
            fseamin = 1.0_dp/fseamin

            write(*, sea)
        end subroutine setup_sea

        subroutine truncate_cplcon_sea()
            call apply_truncation(rhcaps)
            call apply_truncation(rhcapi)
            call apply_truncation(cdsea)
            call apply_truncation(cdice)
            call apply_truncation(beta)
        end subroutine truncate_cplcon_sea
end module
