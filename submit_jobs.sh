#!/bin/bash

for i in {0..4}
do
	sbatch --job-name=SR_EN10s11eMEM$i Q.sh ./run10year.sh SR_EN10s11eMEM$i 06$i 1 10
	sbatch --job-name=SR_CONT10s11eMEM$i Q.sh ./run10year.sh SR_CONT10s11eMEM$i 01$i 0 10
done

for i in {0..4}
do
	sbatch --job-name=SR_RPEFALSE_CONT10s11eMEM$i Q.sh ./run10year.sh SR_RPEFALSE_CONT10s11eMEM$i 01$i 0 10RPEFALSE
done

for i in {0..4}
do
	sbatch --job-name=SR_CONT08s11eMEM$i Q.sh ./run10year.sh SR_CONT08s11eMEM$i 01$i 0 08
	sbatch --job-name=SR_CONT09s11eMEM$i Q.sh ./run10year.sh SR_CONT09s11eMEM$i 01$i 0 09
done

for i in {0..4}
do 
	sbatch --job-name=SR_CONT23s11eMEM$i Q.sh ./run10year.sh SR_CONT23s11eMEM$i 01$i 0 23
done

