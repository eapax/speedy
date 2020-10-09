#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file (0 = no restart)
# $4 = experiment name (directory in which to store experiment)

set -e

if (($# != 4)); then
    echo 'Usage: '${0}' resolution, experiment no., restart no., experiment name' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
TMP=${UT}/output/${4}/exp_${2}
OUT=${UT}/output/${4}/exp_${2}
INP1=${UT}/output/exp_${3}
INP2=${UT}/output/${4}/exp_${3}

mkdir -p ${UT}/output/${4}
mkdir -p ${OUT}

# Setup files
executable=${UT}/source/imp.exe
if ((${3} != 0)); then
    echo "Using restart namelist"
    namelist=${UT}/setup/speedy_default.nml
else
    echo "Using fresh-start namelist"
    namelist=${UT}/setup/speedy_norestart.nml
fi
output=${UT}/setup/default_outputs.nml
precisions=${UT}/setup/10sig11exp.nml

# Copy files from basic version directory
mkdir -p ${TMP}
find ${TMP} -type f -delete
find ${TMP} -type l -delete
find ${TMP} -mindepth 1 -type d -delete
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${output}     ${TMP}/output_requests.nml
cp ${precisions} ${TMP}/precisions.nml

# Link restart file if needed
if ((${3} == 003)); then
  cp ${INP1}/*.rst ${TMP}
elif ((${3} != 0)); then
  cp ${INP2}/*.rst ${TMP}
fi

# Link input files
BC=${UT}/data/bc/${1}
SH=${UT}/hflux

cd ${TMP}
ln -s ${BC}/climatology.nc climatology.nc
#ln -s ${BC}/anomalies.nc   anomalies.nc
#ln -s ${BC}/blank_anomalies.nc   anomalies.nc
ln -s ${BC}/elNinoSSTanomaly_static_nocutoff_scalefactor4_no-mean-adjust.nc anomalies.nc
ln -s ${SH}/hflux_speedy_ver41_1979_2008_clim.grd fort.31

# Link netCDF library
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/netcdf/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/../rpe_complex/lib/

time ./imp.exe | tee out.lis

