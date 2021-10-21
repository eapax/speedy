






#nohup time ./run.sh speedyoneV2_L2_10_RN_1month     010 2 SRoff10   1month > output/speedyoneV2-L2_10_RN_1month.out &
#nohup time ./run.sh speedyoneV2_L2_50_RN_1month     010 2 SRoff50   1month > output/speedyoneV2-L2_50_RN_1month.out &


name = throwaway
nohup time ./run.sh throwaway 010 2 SRoff50   24hour > output/throwaway24.out
cat output/throwaway24.out | mail -s "Process completed" tomkimpson@gmail.com 








