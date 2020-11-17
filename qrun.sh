#!/bin/bash

#sbatch --job-name=EN52s11eMEMBER0 ./newqueue.sh EN52s11eMEMBER0 060 1 52

for i in {0..4}
do
	sbatch --job-name=EN10s11eSRMEM$i ./newqueue.sh EN10s11eSRMEM$i 06$i 1 10
	sbatch --job-name=CONT10s11eSRMEM$i ./newqueue.sh CONT10s11eSRMEM$i 01$i 0 10
done
