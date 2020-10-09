#!/bin/bash

# $1 = experiment name (eg ElNino10s11e)
# $2 = experiment no. for restart file (0 = no restart)
# $3 = choice of SST forcing (0 = no forcing, 1 = El Nino forcing (static, scalefactor = 4, no mean adjust))
# $4 = number of significand bits (52 = double, 23 = single, 10 = half)
# run this job using sbatch with the --job-name=myfavjob flag to set the job name (e.g.=experiment name)

# set the number of nodes
#SBATCH --nodes=1

# set the number of processes per node.
#SBATCH --ntasks-per-node=1

# allocate memory limit
#SBATCH --mem=2000M

# set partition name
#SBATCH --partition=priority-predict,shared

# mail alert at start, end, and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=edmund.paxton@physics.ox.ac.uk

./prun.sh $1 $2 $3 $4
