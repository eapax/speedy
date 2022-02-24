

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
#nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2land_1year                     010 2 SRoff52_cpl1              1year4CO2 > output/speedyoneCPL1.out &
#nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2sea_1year                      010 2 SRoff52_cpl2              1year4CO2 > output/speedyoneCPL2.out &
#nohup time ./run.sh speedyoneCOUPLING_L2_50_RNland2atm_1year                     010 2 SRoff52_cpl3              1year4CO2 > output/speedyoneCPL3.out &
#nohup time ./run.sh speedyoneCOUPLING_L2_50_RNsea2atm_1year                      010 2 SRoff52_cpl4              1year4CO2 > output/speedyoneCPL4.out &

#1 year, 4CO2 breakdown run, agcm
#nohup time ./run.sh speedyoneAGCM_L2_50_RNfordate_1year                     010 2 SRoff52_agcm1              1year4CO2 > output/speedyoneagcm1.out &
#nohup time ./run.sh speedyoneAGCM_L2_50_RNinifluxes_1year                   010 2 SRoff52_agcm2              1year4CO2 > output/speedyoneagcm2.out &
#nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop_1year                      010 2 SRoff52_agcm3              1year4CO2 > output/speedyoneagcm3.out &



#10 year, 4CO2 breakdown run, coupling with a further breakdown around atm2sea
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2land_10year                     010 2 SRoff52_cplA1              10year4CO2 > output/speedyoneCPL1.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNatm2sea_10year                      010 2 SRoff52_cplA2              10year4CO2 > output/speedyoneCPL2.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNland2atm_10year                     010 2 SRoff52_cplA3              10year4CO2 > output/speedyoneCPL3.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNsea2atm_10year                      010 2 SRoff52_cplA4              10year4CO2 > output/speedyoneCPL4.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNinterpolateSST_10year                      010 2 SRoff52_cplA5              10year4CO2 > output/speedyoneCP5L.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNinterpolateSEAICE_10year                   010 2 SRoff52_cplA6              10year4CO2 > output/speedyoneCPL6.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNinterpolateANOM_10year                     010 2 SRoff52_cplA7              10year4CO2 > output/speedyoneCPL7.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNinterpolateOCEAN_10year                    010 2 SRoff52_cplA8              10year4CO2 > output/speedyoneCPL8.out &
# nohup time ./run.sh speedyoneCOUPLING_L2_50_RNinterpolateSEA5_10year                     010 2 SRoff52_cplA9              10year4CO2 > output/speedyoneCPL9.out &


#10 year, 4CO2 breakdown run, agcm with further breakdown on stloop and fordate
# nohup time ./run.sh speedyoneAGCM_L2_50_RNfordate_10year                      010 2 SRoff52_agcmA               10year4CO2 > output/speedyoneagcmA.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNfordate1_10year                     010 2 SRoff52_agcmA1              10year4CO2 > output/speedyoneagcmA1.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNfordate2_10year                     010 2 SRoff52_agcmA2              10year4CO2 > output/speedyoneagcmA2.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNfordate3_10year                     010 2 SRoff52_agcmA3              10year4CO2 > output/speedyoneagcmA3.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop_10year                       010 2 SRoff52_agcmB               10year4CO2 > output/speedyoneagcmB.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop1_10year                      010 2 SRoff52_agcmB1              10year4CO2 > output/speedyoneagcmB1.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop2_10year                      010 2 SRoff52_agcmB2              10year4CO2 > output/speedyoneagcmB2.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop3_10year                      010 2 SRoff52_agcmB3              10year4CO2 > output/speedyoneagcmB3.out &
# nohup time ./run.sh speedyoneAGCM_L2_50_RNstloop4_10year                      010 2 SRoff52_agcmB4              10year4CO2 > output/speedyoneagcmB4.out &

#1 year, 4xCO2 RUN, with precipitation 
#nohup time ./run.sh speedyonePRECIP_L2_52_RN_1year            010 2 SRoff52           1year4CO2 > output/speedyonePRECIP-L2_52_RN_1year.out &
#nohup time ./run.sh speedyonePRECIP_L2_23_RN_1year            010 2 SRoff23           1year4CO2 > output/speedyonePRECIP-L2_23_RN_1year.out &
#nohup time ./run.sh speedyonePRECIP_L2_10_RN_1year            010 2 SRoff10           1year4CO2 > output/speedyonePRECIP-L2_10_RN_1year.out &
#nohup time ./run.sh speedyonePRECIP_L2_10_SR_1year            010 2 10                1year4CO2 > output/speedyonePRECIP-L2_10_SR_1year.out &


#10 year ensemble run

# #Control
# for i in {0..4}
# do
# nohup time ./run.sh speedyoneCONTROL_L2_52_RN_m$i 01$i 2 SRoff52 10year4CO2 > output/speedyoneCONTROL_L2_52_RN_m$i.out &
# done


# #Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITOR_L2_52_RN_m$i 01$i 2 SRoff52 10year4CO2 > output/speedyoneCOMPETITOR_L2_52_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_23_RN_m$i 01$i 2 SRoff23 10year4CO2 > output/speedyoneCOMPETITOR_L2_23_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RN_m$i 01$i 2 SRoff10 10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_SR_m$i 01$i 2 10      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_SR_m$i.out &
# done



# #10 year run to get interpolated land fields
# nohup time ./run.sh speedyoneLAND_L2_52_RN_m1 011 2 SRoff52 10year4CO2 > output/speedyoneLAND_L2_52_RN_m1.out &
# nohup time ./run.sh speedyoneLAND_L2_10_RN_m1 011 2 SRoff10 10year4CO2 > output/speedyoneLAND_L2_10_RN_m1.out &
# nohup time ./run.sh speedyoneLAND_L2_10_SR_m1 011 2 10      10year4CO2 > output/speedyoneLAND_L2_10_SR_m1.out &



#10 year ensemble run extras
#Also do a low precision run, but with a breakdown by sector
#Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNinit_m$i    01$i 2 SRoff10_init    10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNinit_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNagcm_m$i    01$i 2 SRoff10_agcm    10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNagcm_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNcoupler_m$i 01$i 2 SRoff10_coupler 10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNcoupler_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNcheck_m$i   01$i 2 SRoff10_check   10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNcheck_m$i.out &
# done

#10 year ensemble run extras extras
#Also do a low precision run, but with a further breakdown by sector to explore the agcm effects
#Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNfordate_m$i       01$i 2 SRoff10_fordate      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNfordate_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNinifluxes_m$i     01$i 2 SRoff10_inifluxes    10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNinifluxes_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop_m$i        01$i 2 SRoff10_stloop       10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstep_m$i          01$i 2 SRoff10_step         10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstep_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop1_m$i       01$i 2 SRoff10_stloop1      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop1_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop2_m$i       01$i 2 SRoff10_stloop2      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop2_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop3_m$i       01$i 2 SRoff10_stloop3      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop3_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop4_m$i       01$i 2 SRoff10_stloop4      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop4_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNstloop5_m$i       01$i 2 SRoff10_stloop5      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNstloop5_m$i.out &




# done


#10 year ensemble run extras
#Do one more check with agcm and cpl at high prec.
#Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNchecklite_m$i   01$i 2 SRoff10_checklite   10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNchecklite_m$i.out &
# done




# #1 year ensemble run with williams filter turned off

# #Control
# # #Control
# for i in {0..4}
# do
# nohup time ./run.sh speedyoneCONTRORA_L2_52_RN_m$i 01$i 2 SRoff52 1year4CO2 > output/speedyoneCONTROLRA_L2_52_RN_m$i.out &
# done


# #Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITORRA_L2_52_RN_m$i           01$i 2 SRoff52           1year4CO2 > output/speedyoneCOMPETITORRA_L2_52_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA_L2_10_RN_m$i           01$i 2 SRoff10           1year4CO2 > output/speedyoneCOMPETITORRA_L2_10_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA_L2_10_SR_m$i           01$i 2 10                1year4CO2 > output/speedyoneCOMPETITORRA_L2_10_SR_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA_L2_10_RNstloop5_m$i    01$i 2 SRoff10_stloop5   1year4CO2 > output/speedyoneCOMPETITORRA_L2_10_RNstloop5_m$i.out &

# done



# #1 year ensemble run with timeint() modifications to explore truncation theory. H1 = hypothesis one

# #Control
# # #Control
# for i in {0..4}
# do
# nohup time ./run.sh speedyoneCONTROLH1_L2_52_RN_m$i 01$i 2 SRoff52 1year4CO2 > output/speedyoneCONTROLH1_L2_52_RN_m$i.out &
# done


# #Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITORH1_L2_52_RN_m$i           01$i 2 SRoff52           1year4CO2 > output/speedyoneCOMPETITORH1_L2_52_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORH1_L2_10_RN_m$i           01$i 2 SRoff10           1year4CO2 > output/speedyoneCOMPETITORH1_L2_10_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORH1_L2_10_SR_m$i           01$i 2 10                1year4CO2 > output/speedyoneCOMPETITORH1_L2_10_SR_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORH1_L2_10_RNstloop5_m$i    01$i 2 SRoff10_stloop5   1year4CO2 > output/speedyoneCOMPETITORH1_L2_10_RNstloop5_m$i.out &

# done


#10 year ensemble run correct deepdive
#Set coupler and timeint to high prec. THEN get solutions and breakdown.
#Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNboosted_m$i       01$i 2 SRoff10boosted      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNboosted_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITOR_L2_10_SRboosted_m$i       01$i 2 10boosted           10year4CO2 > output/speedyoneCOMPETITOR_L2_10_SRboosted_m$i.out &
# #nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNBstloop1_m$i       01$i 2 SRoff10_stloop1      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNBstloop1_m$i.out &
# #nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNBstloop2_m$i       01$i 2 SRoff10_stloop2      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNBstloop2_m$i.out &
# #nohup time ./run.sh speedyoneCOMPETITOR_L2_10_RNBstloop3_m$i       01$i 2 SRoff10_stloop3      10year4CO2 > output/speedyoneCOMPETITOR_L2_10_RNBstloop3_m$i.out &


# done


# #10 year ensemble run with williams filter turned off

#Control
# #Control
# for i in {0..4}
# do
# nohup time ./run.sh speedyoneCONTRORA10y_L2_52_RN_m$i 01$i 2 SRoff52 10year4CO2 > output/speedyoneCONTROLRA_L2_52_RN_m$i.out &
# done


# #Competitors
# for i in {5..9}
# do
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_52_RN_m$i           01$i 2 SRoff52           10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_52_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_10_RN_m$i           01$i 2 SRoff10           10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_10_RN_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_10_SR_m$i           01$i 2 10                10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_10_SR_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_10_RNstloop5_m$i    01$i 2 SRoff10_stloop5   10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_10_RNstloop5_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_10_RNstloop_m$i    01$i 2 SRoff10_stloop   10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_10_RNstloop_m$i.out &
# nohup time ./run.sh speedyoneCOMPETITORRA10y_L2_10_RNcpl_m$i    01$i 2 SRoff10_coupler   10year4CO2 > output/speedyoneCOMPETITORRA10y_L2_10_RNcpl_m$i.out &


# done


#Try a corrected williams filter

#nohup time ./run.sh speedyoneWILLIAMS_L2_52_RN_10y   010 2 SRoff52               10year4CO2 > output/speedyoneWILLIAMS52_10y.out &
#nohup time ./run.sh speedyoneWILLIAMS_L2_23_RN_10y   010 2 SRoff23               10year4CO2 > output/speedyoneWILLIAMS23_10y.out &
#nohup time ./run.sh speedyoneWILLIAMS_L2_10_RN_10y   010 2 SRoff10               10year4CO2 > output/speedyoneWILLIAMS10_10y.out &
#nohup time ./run.sh speedyoneWILLIAMS_L2_10_SR_10y   010 2 10                    10year4CO2 > output/speedyoneWILLIAMS10SR_10y.out &

#Basic testing for RAW
#nohup time ./run.sh THROWAWAY52   010 2 SRoff52               1year4CO2 > output/THROWAWAY52.out &
#






#100 year ensemble run with correct Williams filtering

# #Control
# for i in {0..4}
# do
# nohup time ./run.sh speedyone100yr_L2_52_RN_m$i 01$i 2 SRoff52 100year4CO2 > output/speedyone100yr_L2_52_RN_m$i.out &
# done


#Competitors - 23/10
#nohup time ./run.sh speedyone100yr_L2_10_SR_m6 016 2 10      100year4CO2 > output/speedyone100yr_L2_10_SR_m6.out &
#nohup time ./run.sh speedyone100yr_L2_23_RN_m6 016 2 SRoff23 100year4CO2 > output/speedyone100yr_L2_23_RN_m6.out &


#Restart experiment




#nohup time ./run.sh speedyone100yr_L2_23_RN_m6_rst 23_m6_rst 2 SRoff23 50year4CO2 > output/speedyone50yr_L2_23_RN_m6_rst.out &
nohup time ./run.sh speedyone100yr_L2_10_SR_m6_rst 10_m6_rst 2 10 83year4CO2 > output/speedyone83yr_L2_10_SR_m6_rst.out &





#for i in {5..9}
#do
# nohup time ./run.sh speedyone100yr_L2_52_RN_m$i 01$i 2 SRoff52 100year4CO2 > output/speedyone100yr_L2_52_RN_m$i.out &
#nohup time ./run.sh speedyone100yr_L2_23_RN_m$i 01$i 2 SRoff23 100year4CO2 > output/speedyone100yr_L2_23_RN_m$i.out &
# nohup time ./run.sh speedyone100yr_L2_10_RN_m$i 01$i 2 SRoff10 100year4CO2 > output/speedyone100yr_L2_10_RN_m$i.out &
# nohup time ./run.sh speedyone100yr_L2_10_SR_m$i 01$i 2 10      100year4CO2 > output/speedyone100yr_L2_10_SR_m$i.out &

# nohup time ./run.sh speedyone100yr_L2_10_RNstloop_m$i        01$i 2 SRoff10_stloop                100year4CO2 > output/speedyone100yr_L2_10_RNstloop_m$i.out &
# nohup time ./run.sh speedyone100yr_L2_10_RNstloop5_m$i       01$i 2 SRoff10_stloop5               100year4CO2 > output/speedyone100yr_L2_10_RNstloop5_m$i.out &
# nohup time ./run.sh speedyone100yr_L2_10_RNcoupler_m$i       01$i 2 SRoff10_coupler               100year4CO2 > output/speedyone100yr_L2_10_RNcoupler_m$i.out &
# nohup time ./run.sh speedyone100yr_L2_10_RNcplstloop5_m$i    01$i 2 SRoff10_cpl_stloop5           100year4CO2 > output/speedyone100yr_L2_10_RNcpl_stloop5_m$i.out &

#done




