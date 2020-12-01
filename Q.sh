#!/bin/bash

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

$*
