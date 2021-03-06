#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of processes per node.
#SBATCH --ntasks-per-node=1

# set max wallclock time (d-h:m:s)

# set partition name
#SBATCH --partition=priority-predict,shared

# set name of job
#SBATCH --job-name=speedy12

# mail alert at start, end, and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=edmund.paxton@physics.ox.ac.uk

nstart=201
nfinish=225
spinup=003

for (( i=$nstart; i<=$nfinish; i++ ))
do
    echo "$i"
    if [ "$i" -eq "$nstart" ]; then
	echo "Restart from $spinup"
	./run12.sh t30 $i $spinup
    else
	echo "Restart from  $((i-1))"
	./run12.sh t30 $i $((i-1))
    fi
    echo "Linking file for next year"
    mv output/exp_$i/1983010100.rst output/exp_$i/1983010100_initial.rst     
    cp -f output/exp_$i/1984010100.rst output/exp_$i/1983010100.rst 
done
