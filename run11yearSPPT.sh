#!/bin/bash

# $1 = experiment name (eg ElNino10s11e)
# $2 = experiment no. for restart file (i.e. initial condition, 0 = from rest)
# $3 = choice of SST anomaly (0 = no forcing, 1 = El Nino (static, scalefactor = 4, no mean adjust))
# $4 = number of significand bits (52 = double, 23 = single, 10 = half)

set -e

#ulimit -n 524288

if (($# != 4)); then
	echo 'Usage: '${0}' experiment name, restart no., choice of SST anomaly forcing (0 = none, 1=El Nino), number of significand bits (52 = double, 23 = single, 10 = half)' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd -P`
TMP=${UT}/output/$1
OUT=/network/group/aopp/predict/TIP016_PAXTON_RPSPEEDY
INP=${UT}/initial_conditions/exp_$2
mkdir -p ${TMP}

# Setup files
executable=${UT}/source/imp.exe
echo "Using restart namelist"
namelist=${UT}/setup/speedy11yearSPPT.nml
output=${UT}/setup/default_outputs.nml
precisions=${UT}/setup/${4}sig11exp.nml

# Copy files from basic version directory
mkdir -p ${TMP}
find ${TMP} -type f -delete
find ${TMP} -type l -delete
find ${TMP} -mindepth 1 -type d -delete
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${output}     ${TMP}/output_requests.nml
cp ${precisions} ${TMP}/precisions.nml

# Link restart file (i.e. set initial conditions)
cp ${INP}/*.rst ${TMP}

# Link input files
BC=${UT}/data/bc/t30
SH=${UT}/hflux

cd ${TMP}
ln -s ${BC}/climatology.nc climatology.nc

#set the SST anomaly forcing
if (($3 == 0)); then
	ln -s ${BC}/blank_anomalies.nc   anomalies.nc
elif (($3 == 1)); then
	ln -s ${BC}/elNinoSSTanomaly_static_nocutoff_scalefactor4_no-mean-adjust.nc anomalies.nc
fi

ln -s ${SH}/hflux_speedy_ver41_1979_2008_clim.grd fort.31

# Link netCDF library
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/netcdf/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/../rpe_complex_stochastic/lib/

time ./imp.exe | tee out.lis

echo 'Run completed, will now move output to shared storage before termination'

cd ${UT}
mv ${TMP} ${OUT}
