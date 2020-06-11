#!/bin/bash
for i in `seq 1 7`;
do
	for j in 1.0 1.5 2.1 3.1 4.6 6.8 9.9
	do
		temp=$((10**$i))
		f_temp=${temp}.0
		val=$(echo "scale=4; $f_temp*$j" | bc)
		echo "Starting to process neutrino energy: $val"
#		cd ../event_extract/
#		python LheToH5.py -e $val
#		cd ../diff_dist/
		python opening_angle.py -e $val &
		echo "Finished processing neutrino energy: $val"
	done
done
wait
echo "All Done!"
