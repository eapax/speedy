#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file (0 = no restart)

set -e

if [ $# -ne 3 ] ; then
    echo 'Usage: '${0}' resolution, experiment no., restart no.' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
TMP=/home/saffin/temp/
OUT=${UT}/output/exp_${2}
INP=${UT}/output/exp_${3}

# Setup files
executable=${UT}/source/imp.exe
namelist=${UT}/setup/speedy_default.nml
output=${UT}/setup/default_outputs.nml
precisions=${UT}/setup/precisions_default.nml

# Copy files from basic version directory
mkdir -p ${TMP}
find ${TMP}/* -type f -delete
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${output}     ${TMP}/output_requests.nml
cp ${precisions} ${TMP}/precisions.nml

# Link restart file if needed
if [ ${3} != 0 ] ; then
  cp ${INP}/*.rst ${TMP}
fi

# Link input files
SB=${UT}/data/bc/${1}/clim
SC=${UT}/data/bc/${1}/anom

cd ${TMP}
ln -sf ${SB}/sfc.grd   fort.20
ln -sf ${SB}/sst.grd   fort.21
ln -sf ${SB}/icec.grd  fort.22
ln -sf ${SB}/stl.grd   fort.23
ln -sf ${SB}/snowd.grd fort.24
ln -sf ${SB}/swet.grd  fort.26
ln -sf ${SC}/ssta.grd  fort.30

ls -l fort.*

# Link rpe shared library
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${UT}/rpe/lib/

time ./imp.exe | tee out.lis

# Copy output to experiment directory
mkdir -p ${OUT}
mv out.lis ${OUT}/atgcm${2}.lis
mv *.rst ${OUT}
mv *.grd ${OUT}
