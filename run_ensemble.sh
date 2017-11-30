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
year='1981'
month='01'
day='01'
hour='00'

# Define directory names
UT=`pwd`
SRC=${UT}/source
TMP=${UT}/tmp
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

if [ $3 != 0 ] ; then
  echo "link restart file ${year}${month}${day}${hour}"
  ln -s ${CD}/${year}${month}${day}${hour}.rst
fi

# Link input files
echo 'link input files to fortran units'
ksh inpfiles.s $1

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
    echo 1 > fort.2
fi
echo ${year} >> fort.2
echo ${month} >> fort.2
echo ${day} >> fort.2
echo ${hour} >> fort.2


# Loop over precisions being tested
for i in {3..23}
do
    echo ${i}
    # Write precision to input file
    # Reduced precision
    echo ${i} > precision.txt
    # Zeroth mode precision
    echo ${i} >> precision.txt
    # Grid-point dynamics precision
    echo ${i} >> precision.txt
    # Initial condition precision
    echo ${i} >> precision.txt

    # Run the model
    time ./imp.exe | tee out.lis

    # Convert .grd output to .nc
    grd2nc.sh yyyymmddhh_p.ctl

    # Move to a unique file labelled by the precision
    mv yyyymmddhh_p.nc ${OUT}/yyyymmddhh_p${i}.nc

    # Remove model output
    mv out.lis ${OUT}/atgcm$2.lis
    rm *.grd
    rm precision.txt
done
