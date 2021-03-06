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
#for i in {1..4}
#do
#	#nohup ./run10year2.sh EN52s11eMEM$i 06$i 1 52 > EN52m$i.out &
#	nohup ./run10year2.sh EN10s11eMEM$i 06$i 1 10 > EN10m$i.out &
#done

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

#nohup time ./run.sh speedtest10 010 2 SRoff10 1year > speedtest10.out & 

#for i in {0..4}
#do
##runonQ	sbatch --job-name=abruptCO2_52sbit_m$i Q.sh ./run_noscratch.sh abrupt4xCO2_52sbit_m$i 01$i 2 52 100year_21ablco2
##runonatm5	nohup time ./run.sh abrupt4xCO2_52sbit_lowCO2_m$i 01$i 2 52 100year > abrupt4xCO2_52sbit_m${i}_lowCO2 &
##runonatm5	nohup time ./run.sh abrupt4xCO2_10sbit_m$i 01$i 2 SRoff10 100year_21ablco2 > abrupt4xCO2_10sbit_m$i &
##nohup time ./run.sh abrupt4xCO2_10sbit_m$i 01$i 2 SRoff10 100year_21ablco2 > abrupt4xCO2_10sbit_m$i.out &
##nohup time ./run.sh abrupt4xCO2_12sbit_m$i 01$i 2 SRoff12 100year_21ablco2 > abrupt4xCO2_12sbit_m$i.out &
#nohup time ./run.sh abrupt4xCO2_10sbitSR_m$i 01$i 2 10 100year_21ablco2 > abrupt4xCO2_10sbitSR_m$i.out &
##runonatm6	nohup time ./run.sh abrupt4xCO2_14sbit_m$i 01$i 2 SRoff14 100year_21ablco2 > abrupt4xCO2_14sbit_m$i &
##runonatm6	nohup time ./run.sh abrupt4xCO2_23sbit_m$i 01$i 2 SRoff23 100year_21ablco2 > abrupt4xCO2_23sbit_m$i &
#done

#rerunning 14 & 23 sbit climate runs (due to lost data on scratch)
#for i in {0..4}
#do
#nohup time ./run.sh abrupt4xCO2_14sbit_m$i 01$i 2 14 65year_21ablco2 > output/abrupt4xCO2_14sbit_m$i.out &
#nohup time ./run.sh abrupt4xCO2_23sbit_m$i 01$i 2 23 65year_21ablco2 > output/abrupt4xCO2_23sbit_m$i.out &
#done

#test that increasing ablco2 indeed increases temperature
#for i in {0..2}
#do
#nohup time ./run.sh test_ablco2_6_m$i 01$i 0 52 1year > output/test_ablco2_6_m$i.out &
#nohup time ./run.sh test_ablco2_10_m$i 01$i 0 52 1year_10ablco2 > output/test_ablco2_10_m$i.out &
#nohup time ./run.sh test_ablco2_20_m$i 01$i 0 52 1year_20ablco2 > output/test_ablco2_10_m$i.out &
#done

#nohup time ./run.sh abrupt4xCO2_10sbitSR_m0_continued 800 2 10 65year_21ablco2_2004 > abrupt4xCO2_10sbitSR_m0_continued.out &
#nohup time ./run.sh abrupt4xCO2_10sbitSR_m1_continued 801 2 10 65year_21ablco2_2004 > abrupt4xCO2_10sbitSR_m1_continued.out &
#nohup time ./run.sh abrupt4xCO2_10sbitSR_m2_continued 802 2 10 65year_21ablco2_2003 > abrupt4xCO2_10sbitSR_m2_continued.out &
#nohup time ./run.sh abrupt4xCO2_10sbitSR_m3_continued 803 2 10 65year_21ablco2_2004 > abrupt4xCO2_10sbitSR_m3_continued.out &
#nohup time ./run.sh abrupt4xCO2_10sbitSR_m4_continued 804 2 10 65year_21ablco2_2006 > abrupt4xCO2_10sbitSR_m4_continued.out &

for i in {0..4}
do 
	nohup time ./run.sh 10sbitSR_m$i 01$i 2 10 65year_21ablco2 > 10sbitSR_m$i.out &
done
