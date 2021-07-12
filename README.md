# speedy
A branch of SPEEDY ver 41, originally cloned from Leo Saffin's gitlab [1], 
which has been ammended to accomodate the reduced precision emulator (RPE) of Dawson & Düben [2]. 

SPEEDY (Simplified Parameterisations primitivE-Equations DYnamics) is a coarse resolution atmospheric circulation model
written by Fred Kucharski & Franco Molteni [3], and this repository has been made public with the permission of Franco Molteni. 
The code here is as cited in the paper [4]. 

## Usage.

To run the model, first edit the `run.sh` script so that `LD_LIBRARY_PATH` points to the RPE library. Next set 
`OUT` to the directory in which you want the output to be 
stored after completion, and `TMP` to the directory in which you want the model to run. 
It's fine to choose the same directory for `OUT` and `TMP`, however for fast performance (especially when using the RPE) 
it's better to set `TMP` to a folder on your computer's scratch space, if possible. 
You can then run

```
./run.sh my_first_experiment 010 0 23 10day
```
The arguments passed to `run.sh` above are as follows:

1: A folder will be created with this name in which your experiment will be contained. 

2: Specifies a choice of initial condition, which in this case is read from the "restart file" contained 
in `initial_conditions/exp_010`. 

3: An integer which specifies the choice of SST anomaly (boundary condition), read from the netcdf file in 
`data/bc/t30/`. The current choices are `0` for no anomaly (i.e. a blank file of zeros), `1` for a constant-in-time field 
which crudely simulates and El Nino year (see [4] for details), and `2` for warming SSTs corresponding to an atmosphere 
subject to an instantaneous quadrupling of CO2 (obtained from a coupled model experiment with EC-Earth)

4: The number of significand bits to be specified by the RPE. Choosing `52` (i.e. double precision) will
turn off the emulator. This argument points the run script to an appropriate namelist file in `setup/precisions`, and not all 
choices of sbits will correspond to a file there (though it is easy to add new namelist files as desired). 

**IMPORTANT:** on the stochastic rounding (SR) branch, by default this currently assumes SR to be turned on! 
To run with SR off one should pass `SRoff23` instead of `23` (see namelist files in `setup/precisions`). 

5: Specifies a namelist file from the `setup` folder prefixed by 'speedy' which contains a list of parameters for the model. 
Primarily this can be used for changing the integration length (for example, `10day` specifies a 10 day integration, with other 
parameters set to some default values), but it can also be used to set more custom namelists (for example, changing the optical thickness of 
atmospheric CO2, or turning on stochastically perturbed parameterisation tendencies (SPPT) as implemented by Saffin [5]). 

## References:

[1] https://gitlab.physics.ox.ac.uk/saffin/speedy

[2] A. Dawson & P. Düben --- _rpe v5: an emulator for reduced floating-point precision in large numerical simulations_ (2017) 

[3] https://www.ictp.it/research/esp/models/speedy.aspx

[4] Paxton et al. --- _Climate Modelling in Low-Precision: Effects of both Deterministic & Stochastic Rounding_ (2021)

[5] Saffin et al. --- _Reduced-precision parametrization: lessons from an intermediate-complexity atmospheric model._ (2020)
