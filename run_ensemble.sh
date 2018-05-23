#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file ( 0 = no restart )


if [ $# -ne 3 ] ; then
    echo 'Usage: '${0}' resol. exp_no. restart_no' 1>&2
    exit 1
fi

# Define precisions to test
pmin=5
pmax=23
# Define directory names
UT=`pwd`
TMP=/home/saffin/temp/
mkdir -p ${UT}/output/exp_${2}
OUT=${UT}/output/exp_${2}
INP=${UT}/output/exp_${3}
CTL=${UT}/output/exp_000

# Setup files
executable=${UT}/source/imp.exe
namelist=${UT}/setup/speedy.nml

# Copy files from basic version directory
mkdir -p ${TMP}
cd ${TMP}
rm *
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${CTL}/*.ctl ${TMP}
cp ${UT}/setup/output_requests.nml ${TMP}
cp ${UT}/setup/precisions.nml ${TMP}

# Link restart file if needed
if [ ${3} != 0 ] ; then
  ln -s ${INP}/*.rst ${TMP}
fi

# Link input files
SB=${UT}/data/bc/${1}/clim
SC=${UT}/data/bc/${1}/anom

ln -sf ${SB}/sfc.grd   fort.20
ln -sf ${SB}/sst.grd   fort.21
ln -sf ${SB}/icec.grd  fort.22
ln -sf ${SB}/stl.grd   fort.23
ln -sf ${SB}/snowd.grd fort.24
ln -sf ${SB}/swet.grd  fort.26
ln -sf ${SC}/ssta.grd  fort.30


# Loop over precisions being tested
for i in $(seq ${pmin} ${pmax})
do
    echo ${i}
    # Write precision to input file
    # reduced_precision
    sed -i "s/rp_convection=.*/rp_convection=${i},/" precisions.nml
    cat precisions.nml

    # Run the model
    time ./imp.exe | tee out.lis

    # Convert .grd output to .nc
    grd2nc_p.sh prognostics_pressure.ctl
    grd2nc.sh tendencies.ctl

    # Move to a unique file labelled by the precision
    mv prognostics_pressure.nc ${OUT}/prognostics_pressure_${i}.nc
    mv tendencies.nc ${OUT}/tendencies_${i}.nc

    # Remove model output
    rm *.grd
done
