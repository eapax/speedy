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

#for i in {0..9}
#do
#	sbatch --job-name=C52MEM$i Q.sh ./run10year.sh CONT52s11eMEM$i 01$i 0 52
#done
#
#for i in {0..4}
#do
#	sbatch --job-name=C50MEM$i Q.sh ./run10year.sh CONT50s11eMEM$i 01$i 0 SRoff50
#done
#

#sbatch --job-name=EN52MEM0 Q.sh ./run10year.sh EN52s11eMEM0 060 1 52

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

#for i in {0..4} 
#do 
#	sbatch --job-name=11m$i Q.sh ./run10year.sh CONT11s11eMEM$i 01$i 0 SRoff11
#	sbatch --job-name=11mSR$i Q.sh ./run10year.sh SR_CONT11s11eMEM$i 01$i 0 11
#	sbatch --job-name=11mSPPT$i Q.sh ./run11yearSPPT.sh SPPT_CONT11s11eMEM$i 01$i 0 SRoff11
#	sbatch --job-name=11mSRSPPT$i Q.sh ./run11yearSPPT.sh SR_SPPT_CONT11s11eMEM$i 01$i 0 11
#done

#for i in {0..3}
#do 
#	sbatch --job-name=10_${i}testspeedy Q.sh ./run10year.sh testCONT10s11eMEM$i 01$i 0 SRoff10 
#done

#sbatch --job-name=EN52speedy Q.sh ./run10year
#sbatch --job-name=EN10speedy Q.sh ./run10year.sh EN10s11eMEM0 060 1 SRoff10

#for i in {0..4}
#do
#	sbatch --job-name=10_${i}speedy Q.sh ./run10year.sh CONT10s11eMEM$i 01$i 0 SRoff10
#	sbatch --job-name=12_${i}speedy Q.sh ./run10year.sh CONT12s11eMEM$i 01$i 0 SRoff12
#	sbatch --job-name=13_${i}speedy Q.sh ./run10year.sh CONT13s11eMEM$i 01$i 0 SRoff13
#	sbatch --job-name=14_${i}speedy Q.sh ./run10year.sh CONT14s11eMEM$i 01$i 0 SRoff14
##	sbatch --job-name=SR14_${i}speedy Q.sh ./run10year.sh SR_CONT14s11eMEM$i 01$i 0 14
#	sbatch --job-name=23_${i}speedy Q.sh ./run10year.sh CONT23s11eMEM$i 01$i 0 SRoff23
#done

#for i in {0..4}
#do
#	sbatch --job-name=alloff10_${i}speedy Q.sh ./run10year.sh alloffCONT10s11eMEM$i 01$i 0 SR+ALLoff10
#	sbatch --job-name=conoff10_${i}speedy Q.sh ./run10year.sh conoffCONT10s11eMEM$i 01$i 0 SR+CONoff10
#done
#
#for i in {0..4}
#do
#	nohup ./run10year2.sh SR_CONT10s11eMEM$i 01$i 0 10 > SR10MEM$i.out &
#	nohup ./run10year2.sh SR_CONT12s11eMEM$i 01$i 0 12 > SR12MEM$i.out &
#	nohup ./run10year2.sh SR_CONT14s11eMEM$i 01$i 0 14 > SR14MEM$i.out &
#done   

#nohup ./run10year2.sh EN52s11eMEM0 060 1 52 > EN52MEM0.out &
#nohup ./run10year2.sh EN10s11eMEM0 060 1 SRoff10 > EN10MEM0.out &

#for i in {0..9}
#do
#	nohup ./run.sh tenday52mem$i 01$i 0 52 10day > tenday52mem$i.out &
#	nohup ./run.sh tenday23mem$i 01$i 0 23 10day > tenday23mem$i.out &
#	nohup ./run.sh tenday52SPPTmem$i 01$i 0 52 10daySPPT > tenday52SPPTmem$i.out &
#done
#
#for i in {0..9}
#do
#	nohup ./run.sh SPPT_CONT52s11eMEM$i 01$i 0 52 11yearSPPT > 11year52SPPTmem$i.out &
#done
#
for i in {1..4}
do
	#nohup ./run10year2.sh EN52s11eMEM$i 06$i 1 52 > EN52m$i.out &
	nohup ./run10year2.sh EN10s11eMEM$i 06$i 1 10 > EN10m$i.out &
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
