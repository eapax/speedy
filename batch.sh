#nohup time ./run.sh training_wheels_23sbit_m0 010 2 SRoff23 10year > output/training_wheels_23sbit_m0.out &
#nohup time ./run.sh training_wheels_14sbit_m0 010 2 SRoff14 10year > output/training_wheels_14sbit_m0.out &


#nohup time ./run.sh training_wheels_12sbit_m0 010 2 SRoff12 10year > output/training_wheels_12sbit_m0.out &
#nohup time ./run.sh training_wheels_10sbit_m0 010 2 SRoff10 10year > output/training_wheels_10sbit_m0.out &
#nohup time ./run.sh training_wheels_10sbit_m0_SR 010 2 10 10year > output/training_wheels_10sbit_m0_SR.out &



#nohup time ./run.sh training_wheels_10sbit_m0 010 2 SRoff10 1year > output/training_wheels_10sbit_m0_1year.out &




#nohup time ./run.sh training_wheels_23sbit_m0_10day 007 2 SRoff23 10day > output/training_wheels_23sbit_m0_10day.out &
#nohup time ./run.sh training_wheels_10sbit_m0_10day 007 2 SRoff10 10day > output/training_wheels_10sbit_m0_10day.out &
#nohup time ./run.sh training_wheels_10sbit_m0_10day_SR 007 2 10 10day > output/training_wheels_10sbit_m0_10day_SR.out &


#El nino
#nohup time ./run.sh el_nino_23sbit_m0_10day 010 1 SRoff23 10day > output/el_nino_23sbit_m0_10day.out &
#nohup time ./run.sh el_nino_wheels_10sbit_m0_10day 010 1 SRoff10 10day > output/el_nino_10sbit_m0_10day.out &
#nohup time ./run.sh el_nino_wheels_10sbit_m0_10day_SR 010 1 10 10day > output/el_nino_10sbit_m0_10day_SR.out &


#El nino static SST. 10 year. ablco2=6
#nohup time ./run.sh el_nino_23sbit_1year 010 1 SRoff23 1year > output/el_nino_23sbit_1year.out &
#nohup time ./run.sh el_nino_10sbit_1year 010 1 SRoff10 1year > output/el_nino_10sbit_1year.out &
#nohup time ./run.sh el_nino_10sbitSR_1year 010 1 10 1year > output/el_nino_10sbit_SR_1year.out &






#ExptCompare - 10days

#L1 - EL NINO BASIC
#nohup time ./run.sh exptcompare_L1_23_RN_dot 010 1 SRoff23 10day > output/exptcompare-L1_23_RN.out &
#nohup time ./run.sh exptcompare_L1_10_RN_dot 010 1 SRoff10 10day > output/exptcompare-L1_10_RN.out &
#nohup time ./run.sh exptcompare_L1_10_SR_dot 010 1 10      10day > output/exptcompare-L1_10_SR.out &

#L2 - L1 + different SST (1 ---> 2)
#nohup time ./run.sh exptcompare_L2_23_RN_dot 010 2 SRoff23 10day > output/exptcompare-L2_23_RN.out &
#nohup time ./run.sh exptcompare_L2_10_RN_dot 010 2 SRoff10 10day > output/exptcompare-L2_10_RN.out &
#nohup time ./run.sh exptcompare_L2_10_SR_dot 010 2 10      10day > output/exptcompare-L2_10_SR.out &


#L3 - L1 + ablco2 level 6---> 21
#nohup time ./run.sh exptcompare_L3_23_RN_dot 010 1 SRoff23 10day_ablco2 > output/exptcompare-L3_23_RN.out &
#nohup time ./run.sh exptcompare_L3_10_RN_dot 010 1 SRoff10 10day_ablco2 > output/exptcompare-L3_10_RN.out &
#nohup time ./run.sh exptcompare_L3_10_SR_dot 010 1 10      10day_ablco2 > output/exptcompare-L3_10_SR.out &


#L4- Both i.e. 4xCO2 run
#nohup time ./run.sh exptcompare_L4_23_RN_dot 010 2 SRoff23 10day_ablco2 > output/exptcompare-L4_23_RN.out &
#nohup time ./run.sh exptcompare_L4_10_RN_dot 010 2 SRoff10 10day_ablco2 > output/exptcompare-L4_10_RN.out &
#nohup time ./run.sh exptcompare_L4_10_SR_dot 010 2 10      10day_ablco2 > output/exptcompare-L4_10_SR.out &


#ExptCompare - 1year. Just comparing L1 and L2

#L1 - EL NINO BASIC
#nohup time ./run.sh exptcompare_L1_23_RN_1year 010 1 SRoff23 1year > output/exptcompare-L1_23_RN.out &
#nohup time ./run.sh exptcompare_L1_10_RN_1year 010 1 SRoff10 1year > output/exptcompare-L1_10_RN.out &
#nohup time ./run.sh exptcompare_L1_10_SR_1year 010 1 10      1year > output/exptcompare-L1_10_SR.out &

#L2 - L1 + different SST (1 ---> 2)
#nohup time ./run.sh exptcompare_L2_23_RN_1year 010 2 SRoff23 1year > output/exptcompare-L2_23_RN_1year.out &
#nohup time ./run.sh exptcompare_L2_10_RN_1year 010 2 SRoff10 1year > output/exptcompare-L2_10_RN_1year.out &
#nohup time ./run.sh exptcompare_L2_10_SR_1year 010 2 10      1year > output/exptcompare-L2_10_SR_1year.out &



#Test I/O with derivatives
#nohup time ./run.sh expt_IO_testing 010 1 SRoff23 10step > output/IO_derivs_test23.out &

#this one with interpolated output levels
#nohup time ./run.sh expt_IO_testing2 010 1 SRoff10 10day > output/IO_derivs102.out &


#this one without interpolated output levels
#nohup time ./run.sh expt_IO_testing3 010 1 SRoff10 10day > output/IO_derivs103.out & 




#Plevels Interpolation Explore - 10days

#L1 - EL NINO BASIC
#nohup time ./run.sh summary_L1_23_RN_dot 010 1 SRoff23 10day > output/summary-L1_23_RN.out &
#nohup time ./run.sh summary_L1_10_RN_dot 010 1 SRoff10 10day > output/summary-L1_10_RN.out &
#nohup time ./run.sh summary_L1_10_SR_dot 010 1 10      10day > output/summary-L1_10_SR.out &

#L2 - L1 + different SST (1 ---> 2)
#nohup time ./run.sh summary_L2_23_RN_dot 010 2 SRoff23 10day > output/summary-L2_23_RN.out &
#nohup time ./run.sh summary_L2_10_RN_dot 010 2 SRoff10 10day > output/summary-L2_10_RN.out &
#nohup time ./run.sh summary_L2_10_SR_dot 010 2 10      10day > output/summary-L2_10_SR.out &


#L3 - L1 + ablco2 level 6---> 21
#nohup time ./run.sh summary_L3_23_RN_dot 010 1 SRoff23 10day_ablco2 > output/summary-L3_23_RN.out &
#nohup time ./run.sh summary_L3_10_RN_dot 010 1 SRoff10 10day_ablco2 > output/summary-L3_10_RN.out &
#nohup time ./run.sh summary_L3_10_SR_dot 010 1 10      10day_ablco2 > output/summary-L3_10_SR.out &


#L4- Both i.e. 4xCO2 run
#nohup time ./run.sh summary_L4_23_RN_dot 010 2 SRoff23 10day_ablco2 > output/summary-L4_23_RN.out &
#nohup time ./run.sh summary_L4_10_RN_dot 010 2 SRoff10 10day_ablco2 > output/summary-L4_10_RN.out &
#nohup time ./run.sh summary_L4_10_SR_dot 010 2 10      10day_ablco2 > output/summary-L4_10_SR.out &




#L1 - EL NINO BASIC
#nohup time ./run.sh lummary_L1_23_RN_dot 010 1 SRoff23 10day > output/lummary-L1_23_RN.out &
#nohup time ./run.sh lummary_L1_10_RN_dot 010 1 SRoff10 10day > output/lummary-L1_10_RN.out &
#nohup time ./run.sh lummary_L1_10_SR_dot 010 1 10      10day > output/lummary-L1_10_SR.out &


#nohup time ./run.sh lummary_L2_23_RN_dot 010 2 SRoff23 10day > output/lummary-L2_23_RN.out &
#nohup time ./run.sh lummary_L2_10_RN_dot 010 2 SRoff10 10day > output/lummary-L2_10_RN.out &
#nohup time ./run.sh lummary_L2_10_SR_dot 010 2 10      10day > output/lummary-L2_10_SR.out &


#L3 - L1 + ablco2 level 6---> 21
#nohup time ./run.sh lummary_L3_23_RN_dot 010 1 SRoff23 10day_ablco2 > output/lummary-L3_23_RN.out &
#nohup time ./run.sh lummary_L3_10_RN_dot 010 1 SRoff10 10day_ablco2 > output/lummary-L3_10_RN.out &
#nohup time ./run.sh lummary_L3_10_SR_dot 010 1 10      10day_ablco2 > output/lummary-L3_10_SR.out &


#L4- Both i.e. 4xCO2 run
#nohup time ./run.sh lummary_L4_23_RN_dot 010 2 SRoff23 10day_ablco2 > output/lummary-L4_23_RN.out &
#nohup time ./run.sh lummary_L4_10_RN_dot 010 2 SRoff10 10day_ablco2 > output/lummary-L4_10_RN.out &
#nohup time ./run.sh lummary_L4_10_SR_dot 010 2 10      10day_ablco2 > output/lummary-L4_10_SR.out &


#L1 - EL NINO BASIC
nohup time ./run.sh geop_L1_50_RN_1year 010 1 SRoff50 1year > output/geop-L1_50_RN.out &
nohup time ./run.sh geop_L1_23_RN_1year 010 1 SRoff23 1year > output/geop-L1_23_RN.out &
nohup time ./run.sh geop_L1_10_RN_1year 010 1 SRoff10 1year > output/geop-L1_10_RN.out &
nohup time ./run.sh geop_L1_10_SR_1year 010 1 10      1year > output/geop-L1_10_SR.out &


nohup time ./run.sh geop_L2_50_RN_1year 010 2 SRoff50 1year > output/geop-L2_50_RN.out &
nohup time ./run.sh geop_L2_23_RN_1year 010 2 SRoff23 1year > output/geop-L2_23_RN.out &
nohup time ./run.sh geop_L2_10_RN_1year 010 2 SRoff10 1year > output/geop-L2_10_RN.out &
nohup time ./run.sh geop_L2_10_SR_1year 010 2 10      1year > output/geop-L2_10_SR.out &


#L3 - L1 + ablco2 level 6---> 21
nohup time ./run.sh geop_L3_50_RN_1year 010 1 SRoff50 1year_ablco2 > output/geop-L3_50_RN.out &
nohup time ./run.sh geop_L3_23_RN_1year 010 1 SRoff23 1year_ablco2 > output/geop-L3_23_RN.out &
nohup time ./run.sh geop_L3_10_RN_1year 010 1 SRoff10 1year_ablco2 > output/geop-L3_10_RN.out &
nohup time ./run.sh geop_L3_10_SR_1year 010 1 10      1year_ablco2 > output/geop-L3_10_SR.out &


#L4- Both i.e. 4xCO2 run
nohup time ./run.sh geop_L4_50_RN_dot 010 2 SRoff50 1year_ablco2 > output/geop-L4_50_RN.out &
nohup time ./run.sh geop_L4_23_RN_dot 010 2 SRoff23 1year_ablco2 > output/geop-L4_23_RN.out &
nohup time ./run.sh geop_L4_10_RN_dot 010 2 SRoff10 1year_ablco2 > output/geop-L4_10_RN.out &
nohup time ./run.sh geop_L4_10_SR_dot 010 2 10      1year_ablco2 > output/geop-L4_10_SR.out &










#test new derivative method
#L1 - EL NINO BASIC
#nohup time ./run.sh expt_tt_sflx_L1_23_RN_dot 010 1 SRoff23 10day > output/expt-L1_23_RN.out &
#nohup time ./run.sh expt_tt_sflx_L1_10_RN_dot 010 1 SRoff10 10day > output/expt-L1_10_RN.out &
#nohup time ./run.sh expt_tt_sflx_L1_10_SR_dot 010 1 10      10day > output/expt-L1_10_SR.out &

#L2 - L1 + different SST (1 ---> 2)
#nohup time ./run.sh expt_tt_sflx_L2_23_RN_dot 010 2 SRoff23 10day > output/expt-L2_23_RN.out &
#nohup time ./run.sh expt_tt_sflx_L2_10_RN_dot 010 2 SRoff10 10day > output/expt-L2_10_RN.out &
#nohup time ./run.sh expt_tt_sflx_L2_10_SR_dot 010 2 10      10day > output/expt-L2_10_SR.out &











