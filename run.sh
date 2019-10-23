#!/usr/bin/env bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file (0 = no restart)

set -e

if ((${#} != 3)); then
    echo 'Usage: '${0}' resolution, experiment no., restart no.' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
TMP=${HOME}/tmp
OUT=${UT}/output/exp_${2}
INP=${UT}/output/exp_${3}

# Setup files
executable=${UT}/source/imp.exe
namelist=${UT}/setup/speedy_default.nml
output=${UT}/setup/default_outputs.nml
precisions=${UT}/setup/precisions_default.nml
start_dump=${INP}/1982010100.rst

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
if ((${3} != 0)); then
  cp -f ${start_dump} ${TMP}
fi

# Link input files
BC=${UT}/data/bc/${1}
SH=${UT}/hflux

cd ${TMP}
ln -s ${BC}/climatology.nc climatology.nc
ln -s ${BC}/anomalies.nc   anomalies.nc
ln -s ${SH}/hflux_speedy_ver41_1979_2008_clim.grd fort.31

# Link netCDF library
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/../rpe_complex/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/netcdf/lib

time ./imp.exe | tee out.lis

# Copy output to experiment directory
mkdir -p ${OUT}
mv *.nc ${OUT}
