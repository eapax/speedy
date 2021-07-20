# speedy
A branch of SPEEDY ver 41, originally cloned from Leo Saffin's gitlab [1], 
which has been ammended to accomodate the reduced precision emulator (RPE) of Dawson & Düben [2]. 
A significant amount of work was done by Leo (last commit in December 2019), 
building on previous work by Sam Hatfield [3], to achieve this, including updating 
the code to a Fortran 90 standard and the introduction of namelist files to simplify the workflow.

SPEEDY (Simplified Parameterisations primitivE-Equations DYnamics) is a coarse resolution atmospheric circulation model
written by Fred Kucharski & Franco Molteni [4] and this repository has been made public with the permission of Franco Molteni. 
The code here is as cited in the paper [5]. Any errors in this code do not reflect those authors.

## Usage.

To run the model, first edit the `run.sh` script so that `LD_LIBRARY_PATH` points to the appropriate RPE library.
Currently, the default branch of this repository is `stochastic_rounding`, which interacts with the branch of the 
RPE developed by Matthew Chantry to support stochastic rounding. 
Next set 
`OUT` to the directory in which you want the output to be 
stored after completion, and `TMP` to the directory in which you want the model to run. 
It's fine to choose the same directory for `OUT` and `TMP`, however for fast performance (especially when using the RPE) 
it's better to set `TMP` to a folder on your computer's scratch space, if possible. 
You can then run the model by

```
./run.sh my_first_experiment 010 0 23 10day
```
The arguments passed to `run.sh` above are as follows (from left to right):

1: A folder will be created with this name in which your experiment will be contained. 

2: Specifies a choice of initial condition, which in this case is read from the "restart file" contained 
in `initial_conditions/exp_010/`. The restart files in `exp_010` through to `exp_019` and in `exp_060`
through to `exp_069` were obtained by
integrating the model for 10 years from rest at 51 significand bits (sbits) of precision (effectively double precision 
plus a tiny rounding error)
with a blank SST anomaly and an El Nino SST anomaly (see below) respectively. 
These can thus be thought of as random initial conditions sampled from the respective attractors.

3: An integer which specifies the choice of SST anomaly (boundary condition), read from the netcdf file in 
`data/bc/t30/`. The current choices are `0` for no anomaly (i.e. a file of zeros), `1` for a constant-in-time field 
which crudely simulates an El Nino year (see [5] for details), and `2` for warming SSTs corresponding to an atmosphere 
subject to an instantaneous quadrupling of CO2 (obtained from a coupled model experiment with EC-Earth). 

4: The number of sbits to be emulated by the RPE. Choosing `52` (i.e. double precision) will
turn off the emulator. Currently the number of exponent bits is 11 (double precision) by default. 
This argument points the run script to an appropriate namelist file in `setup/precisions/`, and not all 
choices of sbits will correspond to a file there (though it is easy to add new namelist files as desired). 

**IMPORTANT:** on the stochastic rounding (SR) branch, by default this currently assumes SR to be turned on! 
To run with SR off one should pass `SRoff23` instead of `23`, see namelist files in `setup/precisions/` for more details. 

A read of the namelist files in `setup/precisions` will show that, in order for the model to be run with 
N sbits across all modules, it is currently necessary to specify the parameter `rp_half_bits=N` as well as 
some additional parameters such as `rp_convection=N`. 
This functionality allows one to selectively turn off low-precision for certain parts of the code, for example in the
convective precipitation parameterisation by setting `rp_convection=52` (or equivalently by omitting to specify
`rp_convection`, which will then default to `52`).  
This functionality is implemented through the `set_precision` subroutine in the `source/mod_prec.f90` module, and 
`set_precision` is currently called only in certain parts of the code, such as convective precipitation.  
One can check which parts of the code it currently affects by running
```
grep -irn set_precision source
```

5: Specifies a namelist file from the `setup/` folder prefixed by `speedy` which contains a list of parameters for the model. 
Primarily this can be used for changing the integration length (for example, `10day` specifies a 10 day integration, with other 
parameters set to some default values), but it can also be used to set more custom namelists (for example, changing the optical thickness of 
atmospheric CO2, or turning on stochastically perturbed parameterisation tendencies (SPPT) as implemented by Saffin [6]). 

## References:

[1] https://gitlab.physics.ox.ac.uk/saffin/speedy

[2] A. Dawson & P. Düben --- _rpe v5: an emulator for reduced floating-point precision in large numerical simulations_ (2017) 

[3] https://samhatfield.co.uk/speedy.f90/

[4] https://www.ictp.it/research/esp/models/speedy.aspx

[5] Paxton et al. --- _Climate Modelling in Low-Precision: Effects of both Deterministic & Stochastic Rounding_ (2021)

[6] Saffin et al. --- _Reduced-precision parametrization: lessons from an intermediate-complexity atmospheric model._ (2020)
