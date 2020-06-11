#!/bin/bash

#Author:Sourav Sarkar
#Objective: This script runs the calchep cross-section calculation
#for different neutrino energies starting from 1GeV upto 99 PeV

#Bash script to run the cross section calculation and loop over if the error is large

#dirname=$1
#mkdir $dirname
#Directory name where to store the cross-section results
dirname=EPA_XS
#parameter to change the number of events
nevents=0
lim1=1.00
lim2=3.00
for i in `seq 0 7`;
do
	for j in 1.0 1.5 2.1 3.1 4.6 6.8 9.9
#	for j in 9.9
	do
		temp=$((10**$i))
                f_temp=${temp}.0
                val=$(echo "scale=4; $f_temp*$j" | bc)
		echo "Working on neutrino energy: $val"
		subdir=En_$val
		mkdir $dirname/$subdir
		q_arr=$(python q_gen.py -e $val)
		for k in $q_arr
		do
			ulim=10.0
			cd results/
		        rm -rf 2_2/
        		echo "Resetting session number..."
        		sed -i "s/#Session_number.*/#Session_number 1/g" session.dat
			echo "Setting photon energy: $k"
			source ../bin/set_momenta $val $k
			echo "Running Initial Vegas Calculation..."
			xsval=$(source ../bin/run_vegas 5 100000 5 100000 &)
			echo $xsval
			xs=$(echo $(cut -d' ' -f1 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
			un=$(echo $(cut -d' ' -f2 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
			count=1
			if [[ $(echo $un'>'$lim1 | bc -l) == "1" ]]; then
				echo "Warning! Error not within limit, trying again..."
			fi
			while [[ $(echo $un'>'$lim1 | bc -l) == "1" && $(echo $un'>'$lim2 | bc -l) == "1" ]];
			do
				echo "Error (%): $un , $ulim"
#			echo "Warning! Error not within limit, trying again..."
				xsval=$(source ../bin/run_vegas 5 100000 5 100000 &)
				xs=$(echo $(cut -d' ' -f1 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
				un=$(echo $(cut -d' ' -f2 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
				if [[ $(echo $un'<'$ulim | bc -l) == "1" ]]; then
					ulim=$un
					echo "Testing $ulim , $un"
				fi
				if [[ "$count" == "6" ]]; then
					lim2=$(echo $ulim'+'1.0 | bc -l)
					echo "Exceeded max #iterations: $count! Trying best available upper limit: $un, $ulim , $lim2"
				fi
				count=$(echo $count'+'1 | bc -l)
			done
			echo "Finished within error tolerance, final error: $un"
#			echo "Starting final calculation..."
#			source ../bin/subproc_cycle $nevents > output.dat &
			echo $val $k $xs $un > data.dat
			wait
#			cp 2_2/run_details.txt ../$dirname/$subdir/xs_$val_$k.txt
			mv data.dat ../$dirname/$subdir/xs_$k.txt
			echo "Moving results to directory complete!"
			echo "Finished calculation for neutrino energy: $val "
			cd ..
done
done
done
