#!/bin/bash
num=0
for i in `seq 5 7`;
do
	for j in 1.0 1.5 2.1 3.1 4.6 6.8 9.9
	do
		temp=$((10**$i))
                f_temp=${temp}.0
                val=$(echo "scale=4; $f_temp*$j" | bc)
		echo "JOB job$num submitter_chep.sub" >> dagman.dag
		echo "CATEGORY job$num GENERIC" >> dagman.dag
		echo "VARS job$num name=\"$num\" args=\"$val\"" >> dagman.dag
		echo "" >> dagman.dag
		num=$(($num+1))
	done
done
