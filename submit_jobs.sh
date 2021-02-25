#!/bin/bash

#for i in {0..4}
#do
#	sbatch --job-name=SR_EN10s11eMEM$i Q.sh ./run10year.sh SR_EN10s11eMEM$i 06$i 1 10
#	sbatch --job-name=SR_CONT10s11eMEM$i Q.sh ./run10year.sh SR_CONT10s11eMEM$i 01$i 0 10
#done
#
#for i in {0..4}
#do
#	sbatch --job-name=SR_RPEFALSE_CONT10s11eMEM$i Q.sh ./run10year.sh SR_RPEFALSE_CONT10s11eMEM$i 01$i 0 10RPEFALSE
#done
#
#for i in {0..4}
#do
#	sbatch --job-name=SR_CONT08s11eMEM$i Q.sh ./run10year.sh SR_CONT08s11eMEM$i 01$i 0 08
#	sbatch --job-name=SR_CONT09s11eMEM$i Q.sh ./run10year.sh SR_CONT09s11eMEM$i 01$i 0 09
#done
#
#for i in {0..4}
#do 
#	sbatch --job-name=SR_CONT23s11eMEM$i Q.sh ./run10year.sh SR_CONT23s11eMEM$i 01$i 0 23
#done
#
#for i in {0..4}
#do
#	nohup ./run10year2.sh CONT12s11eMEM$i 01$i 0 SRoff12 > C12s11eMEM${i}.out &
#	nohup ./run10year2.sh CONT14s11eMEM$i 01$i 0 SRoff14 > C14s11eMEM${i}.out &
#	nohup ./run10year2.sh CONT50s11eMEM$i 01$i 0 SRoff50 > C50s11eMEM${i}.out &
#done

for i in {0..4} 
do 
	sbatch --job-name=11m$i Q.sh ./run10year.sh CONT11s11eMEM$i 01$i 0 SRoff11
	sbatch --job-name=11mSR$i Q.sh ./run10year.sh SR_CONT11s11eMEM$i 01$i 0 11
	sbatch --job-name=11mSPPT$i Q.sh ./run11yearSPPT.sh SPPT_CONT11s11eMEM$i 01$i 0 SRoff11
	sbatch --job-name=11mSRSPPT$i Q.sh ./run11yearSPPT.sh SR_SPPT_CONT11s11eMEM$i 01$i 0 11
done

#for i in {0..9}
#do
#	sbatch --job-name=sppt52s11eMEM$i Q.sh ./run11yearSPPT.sh SPPT_CONT52s11eMEM$i 01$i 0 52
#done
#
#for i in {0..4}
#do
#	sbatch --job-name=sppt10s11eMEM$i Q.sh ./run11yearSPPT.sh SPPT_CONT10s11eMEM$i 01$i 0 SRoff10
#	sbatch --job-name=spptsr10s11eMEM$i Q.sh ./run11yearSPPT.sh SR_SPPT_CONT10s11eMEM$i 01$i 0 10
#	sbatch --job-name=sppt12s11eMEM$i Q.sh ./run11yearSPPT.sh SPPT_CONT12s11eMEM$i 01$i 0 SRoff12
#	sbatch --job-name=spptsr12s11eMEM$i Q.sh ./run11yearSPPT.sh SR_SPPT_CONT12s11eMEM$i 01$i 0 12
#	sbatch --job-name=sppt14s11eMEM$i Q.sh ./run11yearSPPT.sh SPPT_CONT14s11eMEM$i 01$i 0 SRoff14
#	sbatch --job-name=spptsr14s11eMEM$i Q.sh ./run11yearSPPT.sh SR_SPPT_CONT14s11eMEM$i 01$i 0 14 
#	sbatch --job-name=sppt23s11eMEM$i Q.sh ./run11yearSPPT.sh SPPT_CONT23s11eMEM$i 01$i 0 SRoff23
#done

#for i in {0..4}
#do 
#	sbatch --job-name=sr12m$i Q.sh ./run10year.sh SR_CONT12s11eMEM$i 01$i 0 12
#	sbatch --job-name=sr14m$i Q.sh ./run10year.sh SR_CONT14s11eMEM$i 01$i 0 14
#done
