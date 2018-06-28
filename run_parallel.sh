#!/usr/bin/env bash
# $1 = resolution (eg t21, t30)
# $2 = experiment no. for restart file ( 0 = no restart )

if [ $# -ne 2 ] ; then
    echo 'Usage: '$0' resol. restart_no' 1>&2
    exit 1
fi

maxjobs=7
nproc=0

# Precisions
pmin=10
pmax=30

# Start date
day='01'
hour='00'

# Define directory names
UT=`pwd`
SRC=${UT}/source
TMP=${UT}/tmp
CD=${UT}/output/exp_$2

# Copy files from basic version directory
mkdir -p ${TMP}
cd ${TMP}
rm *
cp ${SRC}/yyyymmddhh_p.ctl ${TMP}/
cp ${SRC}/inpfiles.s       ${TMP}/
cp ${SRC}/imp.exe          ${TMP}/

bash inpfiles.s $1

# Loop over start dates
for year in {1987..1991}
do
    for month in `seq -w 01 12`
    do
        cd ${TMP}

        # Link restart dump
        ln -s ${CD}/${year}${month}${day}${hour}.rst

        # Write date input file
        echo 1 > fort.2
        echo ${year} >> fort.2
        echo ${month} >> fort.2
        echo ${day} >> fort.2
        echo ${hour} >> fort.2

        #TDEF 5 LINEAR 00Z01Jan1982 7dy
        DDMMMYYYY=`date -d ${year}-${month}-${day} '+%d%b%Y'`
        echo ${DDMMMYYYY}
        sed -i "8s/.*/TDEF 5 LINEAR 00Z${DDMMMYYYY} 7dy/" yyyymmddhh_p.ctl

        # Loop over precisions being tested
        for i in `seq ${pmin} ${pmax}`
        do
            echo ${i}
            mkdir -p ${TMP}${i}
            cp ${TMP}/* ${TMP}${i}
            cd ${TMP}${i}

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
            ./imp.exe > out.lis 2>&1 &
            nproc=$((nproc+1))

            if [ "$nproc" -ge "$maxjobs" ]; then
                wait
                nproc=0
            fi
        done

        wait
        nproc=0

        for i in `seq ${pmin} ${pmax}`
        do
            echo ${i}

            cd ${TMP}${i}

            # Convert .grd output to .nc
            grd2nc.sh yyyymmddhh_p.ctl

            # Move to a unique file labelled by the precision
            mv yyyymmddhh_p.nc /home/saffin/cirrus/speedy/${year}${month}${day}${hour}_p${i}.nc
            mv out.lis ${UT}/1982-1991/${year}${month}${day}${hour}${i}.lis

            rm *.grd
        done
    done
done


