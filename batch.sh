#L1 - EL NINO BASIC
nohup time ./run.sh speedyone_L1_50_RN_10day 010 1 SRoff50 10day > output/speedyone-L1_50_RN.out &
nohup time ./run.sh speedyone_L1_23_RN_10day 010 1 SRoff23 10day > output/speedyone-L1_23_RN.out &
nohup time ./run.sh speedyone_L1_10_RN_10day 010 1 SRoff10 10day > output/speedyone-L1_10_RN.out &
nohup time ./run.sh speedyone_L1_10_SR_10day 010 1 10      10day > output/speedyone-L1_10_SR.out &


nohup time ./run.sh speedyone_L2_50_RN_10day 010 2 SRoff50 10day > output/speedyone-L2_50_RN.out &
nohup time ./run.sh speedyone_L2_23_RN_10day 010 2 SRoff23 10day > output/speedyone-L2_23_RN.out &
nohup time ./run.sh speedyone_L2_10_RN_10day 010 2 SRoff10 10day > output/speedyone-L2_10_RN.out &
nohup time ./run.sh speedyone_L2_10_SR_10day 010 2 10      10day > output/speedyone-L2_10_SR.out &


#L3 - L1 + ablco2 level 6---> 21
nohup time ./run.sh speedyone_L3_50_RN_10day 010 1 SRoff50 10day_ablco2 > output/speedyone-L3_50_RN.out &
nohup time ./run.sh speedyone_L3_23_RN_10day 010 1 SRoff23 10day_ablco2 > output/speedyone-L3_23_RN.out &
nohup time ./run.sh speedyone_L3_10_RN_10day 010 1 SRoff10 10day_ablco2 > output/speedyone-L3_10_RN.out &
nohup time ./run.sh speedyone_L3_10_SR_10day 010 1 10      10day_ablco2 > output/speedyone-L3_10_SR.out &


#L4- Both i.e. 4xCO2 run
nohup time ./run.sh speedyone_L4_50_RN_10day 010 2 SRoff50 10day_ablco2 > output/speedyone-L4_50_RN.out &
nohup time ./run.sh speedyone_L4_23_RN_10day 010 2 SRoff23 10day_ablco2 > output/speedyone-L4_23_RN.out &
nohup time ./run.sh speedyone_L4_10_RN_10day 010 2 SRoff10 10day_ablco2 > output/speedyone-L4_10_RN.out &
nohup time ./run.sh speedyone_L4_10_SR_10day 010 2 10      10day_ablco2 > output/speedyone-L4_10_SR.out &












