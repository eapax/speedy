#!/bin/bash
nstart=201
nfinish=225
spinup=003
for (( i=$nstart; i<=$nfinish; i++ ))
do
    echo "$i"
    if [ "$i" -eq "$nstart" ]; then
	# cp output/exp_${spinup}/1984010100.rst output/exp_${spinup}/1983010100.rst
	echo "Restart from $spinup"
	./run12.sh t30 $i $spinup
    else
	echo "Restart from  $((i-1))"
	./run12.sh t30 $i $((i-1))
    fi
    echo "Linking file for next year"
    mv output/exp_$i/1983010100.rst output/exp_$i/1983010100_initial.rst     
    cp -f output/exp_$i/1984010100.rst output/exp_$i/1983010100.rst 
done
