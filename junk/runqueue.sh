#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file (0 = no restart)

set -e

if (($# != 3)); then
    echo 'Usage: '${0}' resolution, experiment no., restart no.' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
OUT=${UT}/output/exp_${2}
INP=${UT}/output/exp_${3}

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
precisions=${UT}/setup/precisions_default.nml
queue=${UT}/setup/queue.sh

# Copy files from basic version directory
mkdir -p ${OUT}
find ${OUT} -type f -delete
find ${OUT} -type l -delete
find ${OUT} -mindepth 1 -type d -delete
cp ${executable} ${OUT}/imp.exe
cp ${namelist}   ${OUT}/speedy.nml
cp ${output}     ${OUT}/output_requests.nml
cp ${precisions} ${OUT}/precisions.nml
cp ${queue} ${OUT}/queue.sh


# Link restart file if needed
if ((${3} != 0)); then
  cp ${INP}/*.rst ${OUT}
fi

# Link input files
BC=${UT}/data/bc/${1}
SH=${UT}/hflux

cd ${OUT}
ln -s ${BC}/climatology.nc climatology.nc
#ln -s ${BC}/anomalies.nc   anomalies.nc
#ln -s ${BC}/blank_anomalies.nc   anomalies.nc
ln -s ${BC}/elNino_anomalies.nc anomalies.nc
ln -s ${SH}/hflux_speedy_ver41_1979_2008_clim.grd fort.31

# Link netCDF library
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/netcdf/lib
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/../rpe_complex/lib/


sbatch queue.sh
# time ./imp.exe | tee out.lis

# Copy output to experiment directory
#mkdir -p ${OUT}
#mv model_output.nc ${OUT}

# mv *.rst ${OUT}
