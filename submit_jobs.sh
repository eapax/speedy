#!/bin/bash

#sbatch --job-name=EN52s11eMEMBER0 ./newqueue.sh EN52s11eMEMBER0 060 1 52

#for i in {0..4}
#do
#	sbatch --job-name=SR_EN10s11eMEM$i ./newqueue.sh SR_EN10s11eMEM$i 06$i 1 10
#	sbatch --job-name=SR_CONT10s11eMEM$i ./newqueue.sh SR_CONT10s11eMEM$i 01$i 0 10
#done

for i in {0..4}
do
	sbatch --job-name=SR_RPEFALSE_CONT10s11eMEM$i ./newqueue.sh SR_RPEFALSE_CONT10s11eMEM$i 01$i 0 10RPEFALSE
done

#for i in {0..4}
#do 
#	sbatch --job-name=SR_CONT06s11eMEM$i ./newqueue.sh SR_CONT06s11eMEM$i 01$i 0 06
#done

#for i in {0..4}
#do
#	sbatch --job-name=SR_CONT07s11eMEM$i ./newqueue.sh SR_CONT07s11eMEM$i 01$i 0 07
#done

#for i in {0..4}
#do
#	sbatch --job-name=SR_CONT08s11eMEM$i ./newqueue.sh SR_CONT08s11eMEM$i 01$i 0 08
#	sbatch --job-name=SR_CONT09s11eMEM$i ./newqueue.sh SR_CONT09s11eMEM$i 01$i 0 09
#done

#for i in {0..4}
#do 
#	sbatch --job-name=SR_EN23s11eMEM$i ./newqueue.sh SR_EN23s11eMEM$i 06$i 1 23
#	sbatch --job-name=SR_CONT23s11eMEM$i ./newqueue.sh SR_CONT23s11eMEM$i 01$i 0 23
#done
#
