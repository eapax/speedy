#!/bin/bash

# $1 = experiment name (eg ElNino10s11e)
# $2 = experiment no. for restart file (i.e. initial condition, 0 = from rest)
# $3 = choice of SST anomaly (0 = no forcing, 1 = El Nino (static, scalefactor = 4, no mean adjust), 2 = abrupt 4xCO2 experiment)
# $4 = number of significand bits (52 = double, 23 = single, 10 = half, called from precisions namelist, stochastic rounding on by default unless namelist prefixed with SRoff)
# $5 = namelist (ignore prefix 'speedy' which is implicit) 
#
# example usage: ./run.sh tenday_10sbits_mem$i 01$i 0 SRoff10 10day

set -e

if (($# != 5)); then
	echo 'Usage: '${0}' experiment name, restart no., choice of SST anomaly forcing (0 = none, 1=El Nino, 2 = abrupt 4xCO2 experiment), number of significand bits (eg. 52 = double, 23 = single, 10 = half), namelist (eg. 10daySPPT)' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd -P`
TMP=/local/scratch/kimpson/$1
OUT=/network/group/aopp/predict/TIP016_PAXTON_RPSPEEDY
INP=${UT}/initial_conditions/exp_$2
mkdir -p ${TMP}

# Setup files
executable=${UT}/source/imp.exe
echo "Using restart namelist"
namelist=${UT}/setup/speedy$5.nml
# output=${UT}/setup/default_outputs.nml
#output=${UT}/setup/daily_outputs.nml
output=${UT}/setup/high_frequency_outputs.nml
precisions=${UT}/setup/precisions/${4}sig11exp.nml

# Copy files from basic version directory
find ${TMP} -type f -delete
find ${TMP} -type l -delete
find ${TMP} -mindepth 1 -type d -delete
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${output}     ${TMP}/output_requests.nml
cp ${precisions} ${TMP}/precisions.nml

# Link restart file (i.e. set initial conditions)
if [ $2 != 0 ] ; then
   cp ${INP}/*.rst ${TMP}
fi

# Link input files
BC=${UT}/data/bc/t30
SH=${UT}/hflux

cd ${TMP}
cp ${BC}/climatology.nc climatology.nc

#set the SST anomaly forcing
if (($3 == 0)); then
	cp ${BC}/blank_anomalies.nc anomalies.nc
	echo 'using blank SST anomaly'
elif (($3 == 1)); then
	cp ${BC}/elNinoSSTanomaly_static_nocutoff_scalefactor4_no-mean-adjust.nc anomalies.nc
	echo 'using El Nino SST anomaly'
elif (($3 == 2)); then
	cp ${BC}/abrupt4xCO2anom.nc anomalies.nc
	echo 'using abrupt 4xCO2 SST anomaly'
else
	echo 'usage error: third argument must be 0, 1, or 2. see run script for details.'
fi

cp ${SH}/hflux_speedy_ver41_1979_2008_clim.grd fort.31

# Link netCDF library
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/netcdf/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/../rpe_complex_stochastic/lib/

time ./imp.exe | tee out.lis

echo 'Run completed, will now move output to shared storage before termination'

cd ${UT}
mv ${TMP} ${OUT}
