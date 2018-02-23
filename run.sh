#!/bin/bash

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
TMP=${UT}/tmp
mkdir -p ${UT}/output/exp_${2}
OUT=${UT}/output/exp_${2}
CD=${UT}/output/exp_${3}

# Copy files from basic version directory

echo "copying from ${SRC}/source to ${TMP}"

mkdir -p ${TMP}
cd ${TMP}
rm *

cp ${SRC}/*.f90      ${TMP}/
cp ${SRC}/*.s      ${TMP}/
cp ${SRC}/makefile ${TMP}/
cp ${UT}/setup/speedy.nml ${TMP}/

# Set experiment no. and restart file (if needed)

echo ${3} >  fort.2
echo ${2} >> fort.2

if [ ${3} != 0 ] ; then
  echo "link restart file ${year}${month}${day}${hour}"
  ln -s ${CD}/${year}${month}${day}${hour}.rst
fi 

# Link input files

echo 'link input files to fortran units'

ksh inpfiles.s ${1}

ls -l fort.*

echo ' compiling at_gcm - calling make'

make clean
make imp.exe || { echo "Compilation failed"; exit 1; }

if [ ${4} == make ] ; then
    exit 0
fi

# Write date input file
# First line is 0 for no restart file and 1 for restart
if [ ${3} == 0 ] ; then
    echo 0 > fort.2
else
    echo 1 > fort.2
fi
echo ${year} >> fort.2
echo ${month} >> fort.2
echo ${day} >> fort.2
echo ${hour} >> fort.2

time ./imp.exe | tee out.lis

# Copy output to experiment directory
mv out.lis ${OUT}/atgcm${2}.lis
mv *.rst ${OUT}
mv *.grd ${OUT}
mv *.ctl ${OUT}
