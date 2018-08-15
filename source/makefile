%.o: %.mod

COMPILE=gfortran
COMLIB1=
COMOTT1= -O2
COMCONV= -fconvert=swap

FILES= \
       mod_prec.o \
       cpl_bcinterp.o \
       cpl_land.o \
       cpl_main_interface.o \
       cpl_sea.o \
       cpl_sea_model.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
       dyn_stloop.o \
       ini_agcm_init.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inbcon.o \
       ini_iniatm.o \
       ini_inphys.o \
       ini_invars.o \
       ini_namelist.o \
       ini_rp.o \
       ini_stepone.o \
       interpolation.o \
       mod_atparam.o \
       mod_cli_land.o \
       mod_cli_sea.o \
       mod_cplcon_sea.o \
       mod_cplvar_sea.o \
       mod_cpl_flags.o \
       mod_cpl_land_model.o \
       mod_date.o \
       mod_downscaling.o \
       mod_dyncon0.o \
       mod_dyncon1.o \
       mod_dyncon2.o \
       mod_dynvar.o \
       mod_fft.o \
       mod_fluxes.o \
       mod_fordate.o \
       mod_hdifcon.o \
       mod_physcon.o \
       mod_physvar.o \
       mod_solar.o \
       mod_surfcon.o \
       mod_tsteps.o \
       mod_var_land.o \
       mod_var_sea.o \
       phy_cloud.o \
       phy_convmf.o \
       phy_lscond.o \
       phy_phypar.o \
       phy_radlw.o \
       phy_radsw.o \
       phy_shtorh.o \
       phy_sppt.o \
       phy_suflux.o \
       phy_vdifsc.o \
       ppo_output_stream.o \
       ppo_plevs.o \
       ppo_diagns.o \
       ppo_restart.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o

.SUFFIXES: .f .f90 .o

.f90.o:
	$(COMPILE) $(COMCONV) $(COMOTT1) -c $< -I../rpe/modules

.f.o:
	$(COMPILE) $(COMCONV) $(COMOTT1) -c $<

imp.exe: $(FILES)  at_gcm.o
	$(COMPILE) $(COMCONV) $(COMOTT1) $(FILES) at_gcm.o -o imp.exe $(COMLIB1) -lrpe -L../rpe/lib

.PHONY: clean
clean:
	rm -f *.o *.mod

at_gcm.o               : at_gcm.f90 mod_tsteps.o mod_date.o mod_fordate.o\
                         mod_fluxes.o
cpl_bcinterp.o         : cpl_bcinterp.f90 mod_prec.o
cpl_land.o             : cpl_land.f90 mod_atparam.o mod_var_land.o mod_date.o\
                         mod_cpl_flags.o mod_cpl_land_model.o mod_fluxes.o\
                         mod_cli_land.o mod_var_land.o mod_downscaling.o\
                         mod_prec.o
cpl_main_interface.o   : cpl_main_interface.f90 mod_atparam.o\
                         mod_cpl_land_model.o mod_surfcon.o\
                         mod_cli_land.o mod_cli_sea.o
cpl_sea.o              : cpl_sea.f90 mod_cpl_flags.o mod_atparam.o\
                         mod_cplvar_sea.o mod_date.o mod_fluxes.o\
                         mod_cli_sea.o mod_var_sea.o mod_downscaling.o mod_prec.o
cpl_sea_model.o        : cpl_sea_model.f90 mod_atparam.o mod_cplcon_sea.o\
                         mod_cplvar_sea.o mod_prec.o
dyn_geop.o             : dyn_geop.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o\
                         mod_prec.o
dyn_grtend.o           : dyn_grtend.f90 mod_atparam.o mod_dynvar.o\
                         mod_physvar.o mod_dyncon1.o mod_dyncon2.o\
                         spe_spectral.o mod_prec.o
dyn_implic.o           : dyn_implic.f90 mod_atparam.o mod_dyncon1.o\
                         mod_dyncon2.o mod_prec.o
dyn_sptend.o           : dyn_sptend.f90 mod_atparam.o mod_dynvar.o\
                         mod_dyncon1.o mod_dyncon2.o spe_spectral.o mod_prec.o
dyn_step.o             : dyn_step.f90 mod_dyncon0.o mod_atparam.o\
                         mod_dynvar.o mod_hdifcon.o mod_prec.o
dyn_stloop.o           : dyn_stloop.f90 mod_tsteps.o mod_date.o\
                         mod_dynvar.o phy_radsw.o mod_prec.o
ini_agcm_init.o        : ini_agcm_init.f90 mod_cpl_flags.o mod_date.o\
                         ppo_output_stream.o mod_prec.o mod_fordate.o
ini_impint.o           : ini_impint.f90 mod_dyncon0.o mod_atparam.o\
                         mod_dyncon1.o mod_dyncon2.o mod_hdifcon.o mod_prec.o
ini_inbcon.o           : ini_inbcon.f90 mod_cpl_flags.o mod_date.o\
                         mod_atparam.o mod_surfcon.o mod_cli_land.o\
                         mod_cli_sea.o mod_prec.o spe_spectral.o
ini_indyns.o           : ini_indyns.f90 mod_tsteps.o mod_dyncon0.o\
                         mod_atparam.o mod_dyncon1.o mod_hdifcon.o\
                         spe_spectral.o mod_prec.o
ini_iniatm.o           : ini_iniatm.f90 mod_atparam.o mod_dyncon1.o\
                         mod_fluxes.o mod_prec.o
ini_inphys.o           : ini_inphys.f90 mod_atparam.o mod_physcon.o\
                         phy_convmf.o phy_lscond.o phy_cloud.o phy_radsw.o\
                         phy_radlw.o phy_suflux.o phy_vdifsc.o phy_sppt.o\
                         mod_prec.o
ini_invars.o           : ini_invars.f90 mod_dyncon0.o mod_atparam.o\
                         mod_dynvar.o mod_dyncon1.o mod_surfcon.o mod_date.o\
                         spe_spectral.o mod_prec.o
ini_namelist.o         : ini_namelist.f90 mod_atparam.o mod_downscaling.o\
                         mod_cpl_flags.o mod_cpl_land_model.o mod_cplcon_sea.o\
                         mod_date.o mod_dyncon0.o mod_fordate.o mod_solar.o\
                         mod_surfcon.o mod_tsteps.o phy_cloud.o phy_convmf.o\
                         phy_lscond.o phy_radlw.o phy_radsw.o phy_suflux.o\
                         phy_vdifsc.o phy_sppt.o\
                         mod_cli_land.o mod_cli_sea.o mod_cplvar_sea.o\
                         mod_dyncon1.o mod_dyncon2.o mod_dynvar.o mod_fft.o\
                         mod_fluxes.o mod_hdifcon.o mod_physcon.o mod_physvar.o\
                         mod_var_land.o mod_var_sea.o spe_spectral.o
ini_rp.o               : ini_rp.f90 mod_dyncon1.o mod_physcon.o mod_prec.o\
                         mod_cli_land.o mod_cli_sea.o mod_cpl_land_model.o\
                         mod_cplcon_sea.o mod_cplvar_sea.o mod_dyncon0.o\
                         mod_dyncon1.o mod_dyncon2.o mod_fft.o mod_fordate.o\
                         mod_hdifcon.o mod_physcon.o mod_solar.o mod_surfcon.o\
                         mod_tsteps.o mod_var_land.o mod_var_sea.o phy_convmf.o\
                         phy_lscond.o phy_cloud.o phy_radsw.o phy_radlw.o\
                         phy_suflux.o phy_vdifsc.o phy_sppt.o spe_spectral.o
ini_stepone.o          : ini_stepone.f90 mod_tsteps.o mod_date.o mod_prec.o
interpolation.o        : interpolation.f90 mod_prec.o
mod_atparam.o          : mod_atparam.f90
mod_cli_land.o         : mod_cli_land.f90 mod_atparam.o
mod_cli_sea.o          : mod_cli_sea.f90 mod_atparam.o
mod_cpl_flags.o        : mod_cpl_flags.f90
mod_cpl_land_model.o   : mod_cpl_land_model.f90 mod_atparam.o mod_prec.o
mod_cplcon_sea.o       : mod_cplcon_sea.f90 mod_atparam.o mod_prec.o
mod_cplvar_sea.o       : mod_cplvar_sea.f90 mod_atparam.o
mod_date.o             : mod_date.f90 mod_atparam.o mod_tsteps.o mod_prec.o
mod_downscaling.o      : mod_downscaling.f90 mod_atparam.o spe_spectral.o\
                         interpolation.o mod_prec.o
mod_dyncon0.o          : mod_dyncon0.f90 mod_prec.o
mod_dyncon1.o          : mod_dyncon1.f90 mod_atparam.o mod_prec.o
mod_dyncon2.o          : mod_dyncon2.f90 mod_atparam.o
mod_dynvar.o           : mod_dynvar.f90 mod_atparam.o
mod_fft.o              : mod_fft.f90 mod_atparam.o
mod_fluxes.o           : mod_fluxes.f90 mod_atparam.o mod_prec.o mod_date.o\
                         mod_tsteps.o mod_physcon.o mod_surfcon.o mod_physvar.o\
                         mod_var_sea.o mod_fordate.o phy_radlw.o
mod_fordate.o          : mod_fordate.f90 mod_prec.o mod_dyncon0.o \
                         mod_atparam.o mod_hdifcon.o mod_physcon.o\
                         mod_surfcon.o mod_cli_land.o mod_date.o\
                         mod_var_land.o mod_cli_sea.o mod_var_sea.o\
                         mod_solar.o phy_shtorh.o\
                         spe_spectral.o
mod_hdifcon.o          : mod_hdifcon.f90 mod_atparam.o
mod_physcon.o          : mod_physcon.f90 mod_atparam.o mod_prec.o
mod_physvar.o          : mod_physvar.f90 mod_atparam.o
mod_prec.o             : mod_prec.f90
mod_solar.o            : mod_solar.f90 mod_atparam.o mod_prec.o mod_physcon.o
mod_surfcon.o          : mod_surfcon.f90 mod_atparam.o mod_prec.o
mod_tsteps.o           : mod_tsteps.f90 mod_prec.o
mod_var_land.o         : mod_var_land.f90 mod_atparam.o
mod_var_sea.o          : mod_var_sea.f90 mod_atparam.o
phy_cloud              : phy_cloud.f90 mod_atparam.o mod_prec.o
phy_convmf.o           : phy_convmf.f90 mod_atparam.o mod_physcon.o mod_prec.o
phy_lscond.o           : phy_lscond.f90 mod_atparam.o mod_physcon.o mod_prec.o
phy_phypar.o           : phy_phypar.f90 mod_atparam.o mod_physvar.o\
                         mod_cpl_flags.o mod_physcon.o mod_surfcon.o\
                         mod_var_land.o mod_var_sea.o mod_fluxes.o\
                         phy_shtorh.o phy_convmf.o phy_lscond.o phy_cloud.o\
                         phy_radsw.o phy_radlw.o phy_suflux.o phy_vdifsc.o\
                         phy_sppt.o mod_prec.o
phy_radlw.o            : phy_radlw.f90 mod_atparam.o mod_prec.o mod_physcon.o\
                         mod_physvar.o
phy_radsw.o            : phy_radsw.f90 mod_atparam.o mod_prec.o mod_physcon.o\
                         mod_physvar.o mod_fordate.o mod_solar.o phy_radlw.o
phy_shtorh.o           : phy_shtorh.f90 mod_prec.o
phy_sppt.o             : mod_atparam.o mod_prec.o mod_tsteps.o mod_dyncon1.o\
                         spe_spectral.o
phy_suflux.o           : phy_suflux.f90 mod_atparam.o mod_prec.o mod_physcon.o\
                         mod_fordate.o phy_radlw.o phy_shtorh.o mod_surfcon.o
phy_vdifsc.o           : phy_vdifsc.f90 mod_atparam.o mod_prec.o mod_physcon.o
ppo_diagns.o           : ppo_diagns.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o\
                         ppo_output_stream.o spe_spectral.o mod_prec.o
ppo_output_stream.o    : ppo_output_stream.f90 mod_atparam.o mod_dynvar.o\
                         mod_physvar.o mod_physcon.o mod_date.o ppo_plevs.o\
                         spe_spectral.o mod_prec.o
ppo_plevs.o            : ppo_plevs.f90 mod_atparam.o mod_dynvar.o mod_physvar.o\
                         mod_physcon.o spe_spectral.o mod_prec.o
ppo_restart.o          : ppo_restart.f90 mod_atparam.o mod_dynvar.o mod_date.o\
                         mod_downscaling.o mod_prec.o
spe_matinv.o           : spe_matinv.f90 mod_prec.o
spe_spectral.o         : spe_spectral.f90 mod_atparam.o mod_prec.o
spe_subfft_fftpack.o   : spe_subfft_fftpack.f90 mod_atparam.o spe_spectral.o\
                         mod_fft.o mod_prec.o
