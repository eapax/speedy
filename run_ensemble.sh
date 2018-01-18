#!/usr/bin/env bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file ( 0 = no restart )
# $4 = make only, and don't run? ("make" for yes, "run" for no)


if [ $# -ne 4 ] ; then
    echo 'Usage: '$0' resol. exp_no. restart_no make_only_or_run' 1>&2
    exit 1
fi

# Start date
year='1982'
month='01'
day='01'
hour='00'

# Define directory names
UT=`pwd`
SRC=${UT}/source
TMP=${UT}/temp
mkdir -p ${UT}/output/exp_$2
OUT=${UT}/output/exp_$2
CD=${UT}/output/exp_$3

# Copy files from basic version directory
echo "copying from $SRC/source to $TMP"
mkdir -p ${TMP}
cd ${TMP}
rm *
cp ${OUT}/yyyymmddhh_p.ctl ${TMP}/
cp ${SRC}/*.f90      ${TMP}/
cp ${SRC}/*.h      ${TMP}/
cp ${SRC}/*.s      ${TMP}/
cp ${SRC}/makefile ${TMP}/

# Link input files
echo 'link input files to fortran units'
sh inpfiles.s $1

ls -l fort.*

echo ' compiling at_gcm - calling make'

make clean
make imp.exe || { echo "Compilation failed"; exit 1; }

if [ $4 == make ] ; then
    exit 0
fi

# Write date input file
if [ $3 = 0 ]; then
    echo 0 > fort.2
else
    echo "link restart file ${year}${month}${day}${hour}"
    ln -s ${CD}/${year}${month}${day}${hour}.rst
    echo 1 > fort.2
fi
echo ${year} >> fort.2
echo ${month} >> fort.2
echo ${day} >> fort.2
echo ${hour} >> fort.2


# Loop over precisions being tested
for i in {5..52}
do
    echo ${i}
    # Write precision to input file
    # reduced_precision
    echo ${i} > precision.txt
    # initial values
    echo ${i} >> precision.txt
    # spectral transform
    echo ${i} >> precision.txt
    # grid physics
    echo ${i} >> precision.txt
    # convection
    echo ${i} >> precision.txt
    # condensation
    echo ${i} >> precision.txt
    # short-wave radiation
    echo ${i} >> precision.txt
    # long-wave radiation
    echo ${i} >> precision.txt
    # surface fluxes
    echo ${i} >> precision.txt
    # vertical diffusion
    echo ${i} >> precision.txt
    # SPPT
    echo ${i} >> precision.txt
    # grid dynamics
    echo ${i} >> precision.txt
    # spectral dynamics
    echo ${i} >> precision.txt
    # diffusion
    echo ${i} >> precision.txt
    # time stepping
    echo ${i} >> precision.txt
    # Prognostics
    echo ${i} >> precision.txt
    # Tendencies
    echo ${i} >> precision.txt
    # Initialisation
    echo ${i} >> precision.txt
    # Parameters
    echo ${i} >> precision.txt


    # Run the model
    time ./imp.exe | tee out.lis

    # Convert .grd output to .nc
    grd2nc_p.sh yyyymmddhh_p.ctl

    # Move to a unique file labelled by the precision
    mv yyyymmddhh_p.nc ${OUT}/yyyymmddhh_p${i}.nc

    # Remove model output
    mv out.lis ${OUT}/atgcm$2.lis
    rm *.grd
    rm precision.txt
done
