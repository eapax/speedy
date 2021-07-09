! ug1/vg1 (ms-1)
case(1)
name = 'zonal_velocity'
units = 'm s-1'

! vg1
case(2)
name = 'meridional_velocity'
units = 'm s-1'

! tg1 (K)
case(3)
name = 'temperature'
units = 'K'

! qg1 (kg/kg)
case(4)
name = 'specific_humidity'
units = ''

! phig1  = geopotential
case(5)
name = 'geopotential_height'
units = 'm'

! pslg1
case(6)
name = 'logarithm_of_surface_pressure'
units = ''
l_3d = .false.

! se     = dry static energy
case(7)
name = 'dry_static_energy'
units = 'J kg-1'

! rh     = relative humidity
case(8)
name = 'relative_humidity'
units = ''

! qsat   = saturation specific humidity (g/kg)
case(9)
name = 'saturation_specific_humidity'
units = 'g kg-1'

! psg    = surface pressure
case(10)
name = 'surface_pressure'
units = 'Pa'
l_3d = .false.

! ts     = surface temperature
case(11)
name = 'surface_temperature'
units = 'K'
l_3d = .false.

! tskin  = skin temperature
case(12)
name = 'skin_temperature'
units = 'K'
l_3d = .false.

! u0     = near-surface u-wind
case(13)
name = 'near_surface_zonal_velocity'
units = 'm s-1'
l_3d = .false.

! v0     = near-surface v-wind
case(14)
name = 'near_surface_meridional_velocity'
units = 'm s-1'
l_3d = .false.

! t0     = near-surface air temperature
case(15)
name = 'near_surface_temperature'
units = 'K'
l_3d = .false.

! q0     = near-surface specific humidity (g/kg)
case(16)
name = 'near_surface_specific_humidity'
units = 'g kg-1'
l_3d = .false.

! cloudc = total cloud cover (fraction)
case(17)
name = 'cloud_cover'
units = ''
l_3d = .false.

! clstr  = stratiform cloud cover (fraction)
case(18)
name = 'stratiform_cloud_cover'
units = ''
l_3d = .false.

! cltop  = norm. pressure at cloud top
case(19)
name = 'cloud_top_pressure'
units = ''
l_3d = .false.

! prtop  = top of precipitation (level index)
case(20)
name = 'level_of_precipitation'
units = ''
l_3d = .false.

! precnv = convective precipitation  [g/(m^2 s)], total
case(31)
name = 'convective_precipitation'
units = 'g m-2 s-1'
l_3d = .false.

! precls = large-scale precipitation [g/(m^2 s)], total
case(32)
name = 'large_scale_precipitation'
units = 'g m-2 s-1'
l_3d = .false.

! snowcv = convective precipitation  [g/(m^2 s)], snow only
!case(33)

! snowls = large-scale precipitation [g/(m^2 s)], snow only
!case(34)

! cbmf   = cloud-base mass flux
case(35)
name = 'cloud_base_mass_flux'
units ='unknown'
l_3d = .false.

! tsr    = top-of-atm. shortwave radiation (downward)
case(36)
name = 'top_of_atmosphere_shortwave_radiation'
units = 'unknown'
l_3d = .false.

! ssrd   = surface shortwave radiation (downward-only)
case(37)
name = 'downward_shortwave_radiation_at_surface'
units = 'unknown'
l_3d = .false.

! ssr    = surface shortwave radiation (net downward)
case(38)
name = 'net_downward_shortwave_radiation_at_surface'
units = 'unknown'
l_3d = .false.

! slrd   = surface longwave radiation  (downward-only)
case(39)
name = 'downward_longwave_radiation_at_surface'
units = 'unknown'
l_3d = .false.

! slr    = surface longwave radiation  (net upward)
case(40)
name = 'net_upward_longwave_radiation_at_surface'
units = 'unknown'
l_3d = .false.

! olr    = outgoing longwave radiation (upward)
case(41)
name = 'outgoing_longwave_radiation'
units = 'unknown'
l_3d = .false.

! slru   = surface longwave emission   (upward)
!                                   (1:land, 2:sea, 3: wgt. average)
case(42)
name = 'surface_longwave_emission'
units = 'unknown'
l_3d = .false.

! ustr   = u-stress                 (1:land, 2:sea, 3: wgt. average)
case(43)
name = 'zonal_wind_stress'
units = 'unknown'
l_3d = .false.

! vstr   = v-stress                 (1:land, 2:sea, 3: wgt. average)
case(44)
name = 'meridional_wind_stress'
units = 'unknown'
l_3d = .false.

! shf    = sensible heat flux       (1:land, 2:sea, 3: wgt. average)
case(45)
name = 'sensible_heat_flux'
units = 'unknown'
l_3d = .false.

! evap   = evaporation [g/(m^2 s)]  (1:land, 2:sea, 3: wgt. average)
case(46)
name = 'evaporation'
units = 'g m-2 s-1'
l_3d = .false.

! hfluxn = net heat flux into surf. (1:land, 2:sea, 3: ice-sea dif.)
case(47)
name = 'net_heat_flux_into_surface'
units = ''
l_3d = .false.

! alb_l = Surface albedo over land
case(48)
name = 'surface_albedo_over_land'
units = 'unknown'
l_3d = .false.

! alb_s = surface albedo over sea
case(49)
name = 'surface_albedo_over_sea'
units = 'unknown'
l_3d = .false.

! albsfc = surface albedo
case(50)
name = 'surface_albedo'
units = 'unknown'
l_3d = .false.

! snowc = snow cover
case(51)
name = 'snow_cover'
units = 'unknown'
l_3d = .false.

! fsol = flux of incoming solar radiation
case(52)
name = 'flux_of_incoming_solar_radiation'
units = 'unknown'
l_3d = .false.

! ozone = ozone concentration
case(53)
name = 'ozone'
units = 'unknown'
l_3d = .false.

! ozupp = ozone concentration ...
case(54)
name = 'ozupp'
units = 'unknown'
l_3d = .false.

! zenit = ???
case(55)
name = 'zenit'
units = 'unknown'
l_3d = .false.

! stratz = ???
case(56)
name = 'stratz'
units = 'unknown'
l_3d = .false.

! stl_am = Surface temperature from the land model
case(57)
name = 'land_model_surface_temperature'
units = 'K'
l_3d = .false.

! soilw_am = Soil water from the land model
case(58)
name = 'land_model_soil_moisture'
units = 'unknown'
l_3d = .false.

! sst_am = Sea surface temperature from the ocean model
case(59)
name = 'ocean_model_surface_temperature'
units = 'unknown'
l_3d = .false.

! ssti_om = Sea surface + ice temperature from the ocean model
case(60)
name = 'ocean_and_ice_model_surface_temperature'
units = 'K'
l_3d = .false.

! phis0 = spectrally filtered surface geopotential
case(61)
name = 'surface_geopotential'
units = 'unknown'
l_3d = .false.

! fmask1 = fractional land-sea mask
case(62)
name = 'fractional_land_sea_mask'
units = ''
l_3d = .false.

! tt_cnv  =  temperature tendency due to convection
case(101)
name = 'temperature_tendency_due_to_convection'
units = 'K s-1'

! qt_cnv  = sp. humidity tendency due to convection
case(102)
name = 'specific_humidity_tendency_due_to_convection'
units = 's-1'

! tt_lsc  =  temperature tendency due to large-scale condensation
case(103)
name = 'temperature_tendency_due_to_condensation'
units = 'K s-1'

! qt_lsc  = sp. humidity tendency due to large-scale condensation
case(104)
name = 'specific_humidity_tendency_due_to_condensation'
units = 's-1'

! tt_rsw  =  temperature tendency due to short-wave radiation
case(105)
name = 'temperature_tendency_due_to_shortwave_radiation'
units = 'K s-1'

! tt_rlw  =  temperature tendency due to long-wave radiation
case(106)
name = 'temperature_tendency_due_to_longwave_radiation'
units = 'K s-1'

! ut_sflx  =       u-wind tendency due to surface fluxes
case(107)
name = 'zonal_velocity_tendency_due_to_surface_fluxes'
units = 'm s-2'

! vt_sflx  =       v-wind tendency due to surface fluxes
case(108)
name = 'meridional_velocity_tendency_due_to_surface_fluxes'
units = 'm s-2'

! tt_sflx  =  temperature tendency due to surface fluxes
case(109)
name = 'temperature_tendency_due_to_surface_fluxes'
units = 'K s-1'

! qt_sflx  = sp. humidity tendency due to surface fluxes
case(110)
name = 'specific_humidity_tendency_due_to_surface_fluxes'
units = 's-1'

! ut_pbl  =       u-wind tendency due to PBL and diffusive processes
case(111)
name = 'zonal_velocity_tendency_due_to_vertical_diffusion'
units = 'm s-2'

! vt_pbl  =       v-wind tendency due to PBL and diffusive processes
case(112)
name = 'meridional_velocity_tendency_due_to_vertical_diffusion'
units = 'm s-2'

! tt_pbl  =  temperature tendency due to PBL and diffusive processes
case(113)
name = 'temperature_tendency_due_to_vertical_diffusion'
units = 'K s-1'

! qt_pbl  = sp. humidity tendency due to PBL and diffusive processes
case(114)
name = 'specific_humidity_tendency_due_to_vertical_diffusion'
units = 's-1'

! ut_phy  =       u-wind tendency due to all physics processes
case(115)
name = 'zonal_velocity_tendency_due_to_all_parametrizations'
units = 'm s-2'

! vt_phy  =       v-wind tendency due to all physics processes
case(116)
name = 'meridional_velocity_tendency_due_to_all_parametrizations'
units = 'm s-2'

! tt_phy  =  temperature tendency due to all physics processes
case(117)
name = 'temperature_tendency_due_to_all_parametrizations'
units = 'K s-1'

! qt_phy  = sp. humidity tendency due to all physics processes
case(118)
name = 'specific_humidity_tendency_due_to_all_parametrizations'
units = 's-1'

! ut_sppt =       u-wind tendency due to stochastic perturbation
case(119)
name = 'zonal_velocity_tendency_due_to_stochastic_perturbation'
units = 'm s-2'

! vt_sppt =       v-wind tendency due to stochastic perturbation
case(120)
name = 'meridional_velocity_tendency_due_to_stochastic_perturbation'
units = 'm s-2'

! tt_sppt =  temperature tendency due to stochastic perturbation
case(121)
name = 'temperature_tendency_due_to_stochastic_perturbation'
units = 'K s-1'

! qt_sppt = sp. humidity tendency due to stochastic perturbation
case(122)
name = 'specific_humidity_tendency_due_to_stochastic_perturbation'
units = 's-1'

! 3D Stochastic perturbation pattern
case(123)
name = 'stochastic_perturbation'
units = ''


