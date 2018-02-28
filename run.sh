#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file ( 0 = no restart )


if [ $# -ne 3 ] ; then
    echo 'Usage: '${0}' resol. exp_no. restart_no' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
TMP=/home/saffin/temp/
mkdir -p ${UT}/output/exp_${2}
OUT=${UT}/output/exp_${2}
INP=${UT}/output/exp_${3}

# Start date
year='1982'
month='01'
day='01'
hour='00'

# Setup files
executable=${UT}/source/imp.exe
namelist=${UT}/setup/speedy.nml


# Copy files from basic version directory
mkdir -p ${TMP}
cd ${TMP}
rm *
cp ${executable} ${TMP}/imp.exe
cp ${namelist}   ${TMP}/speedy.nml
cp ${UT}/setup/output_requests.nml ${TMP}

# Link restart file if needed
if [ ${3} != 0 ] ; then
  echo "link restart file ${year}${month}${day}${hour}"
  ln -s ${INP}/${year}${month}${day}${hour}.rst
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

ls -l fort.*

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
