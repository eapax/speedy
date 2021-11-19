

#10 year, 4xCO2 RUN 
#nohup time ./run.sh speedyoneFIG1_L2_52_RN_10year            010 2 SRoff52           10year4CO2 > output/speedyoneFIG1-L2_52_RN_10year.out &
#nohup time ./run.sh speedyoneFIG1_L2_23_RN_10year            010 2 SRoff23           10year4CO2 > output/speedyoneFIG1-L2_23_RN_10year.out &
#nohup time ./run.sh speedyoneFIG1_L2_10_RN_10year            010 2 SRoff10           10year4CO2 > output/speedyoneFIG1-L2_10_RN_10year.out &
#nohup time ./run.sh speedyoneFIG1_L2_52_SR_10year            010 2 52                10year4CO2 > output/speedyoneFIG1-L2_52_SR_10year.out &
#nohup time ./run.sh speedyoneFIG1_L2_23_SR_10year            010 2 23                10year4CO2 > output/speedyoneFIG1-L2_23_SR_10year.out &
#nohup time ./run.sh speedyoneFIG1_L2_10_SR_10year            010 2 10                10year4CO2 > output/speedyoneFIG1-L2_10_SR_10year.out &


#nohup time ./run.sh speedyoneFIG1_L2_10_SR_10year                   010 2 10                10year4CO2 > output/speedyoneFIG1-L2_10_SR_10year.out &


#10 year, 4xCO2 breakdown run, umbrella
#nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNinit_10year                      010 2 SRoff52_init              10year4CO2 > output/speedyoneBREAKDOWN-L2_50_RNinit_10year.out &
#nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNatmosphere_10year                010 2 SRoff52_agcm              10year4CO2 > output/speedyoneBREAKDOWN-L2_50_RNagcm_10year.out &
#nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNcoupling_10year                  010 2 SRoff52_coupler           10year4CO2 > output/speedyoneBREAKDOWN-L2_50_RNcoupling_10year.out &



#1 year, 4CO2 breakdown run, coupling
nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2land_1year                     010 2 SRoff52_cpl1              1year4CO2 > output/speedyoneCPL1.out &
nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2sea_1year                      010 2 SRoff52_cpl2              1year4CO2 > output/speedyoneCPL2.out &
nohup time ./run.sh speedyoneCOUPLING_L2_50_RNland2atm_1year                     010 2 SRoff52_cpl3              1year4CO2 > output/speedyoneCPL3.out &
nohup time ./run.sh speedyoneCOUPLING_L2_50_RNsea2atm_1year                      010 2 SRoff52_cpl4              1year4CO2 > output/speedyoneCPL4.out &




# nohup time ./run.sh speedyoneBREAKDOWN_L2_10_RN_1year                  010 2 SRoff10                   1year > output/speedyoneBREAKDOWN-L2_10_RN_1year.out &
# nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNcoupler_1year           010 2 SRoff52_coupler           1year > output/speedyoneBREAKDOWN-L2_50_RNcoupler_1year.out &
# nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNagcm_1year              010 2 SRoff52_agcm              1year > output/speedyoneBREAKDOWN-L2_50_RNagcm_1year.out &
# nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNinit_1year              010 2 SRoff52_init              1year > output/speedyoneBREAKDOWN-L2_50_RNinit_1year.out &
# nohup time ./run.sh speedyoneBREAKDOWN_L2_50_RNinterpolation_1year     010 2 SRoff52_interpolation     1year > output/speedyoneBREAKDOWN-L2_50_RNinterpolation_1year.out &
# nohup time ./run.sh speedyoneBREAKDOWN_L2_52_RN_1year                  010 2 SRoff52                   1year > output/speedyoneBREAKDOWN-L2_52_RN_1year.out &

