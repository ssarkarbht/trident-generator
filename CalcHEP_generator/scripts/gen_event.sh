#!/bin/sh
#Script to run CalcHEP event generation

target_path=/data/user/ssarkar/CLHP_TEST/test_dir
current_path=$(eval pwd)
dir=$target_path/En_$1
source_dir=/data/user/ssarkar/CalcHEP/CHDF_event/calchep_3.7.1
#parameter to change the number of events
#nevents=100000
nevents=100
lim1=1.00
lim2=3.00
source $source_dir/mkWORKdir $dir
./calchep -blind "[[[[[[[[[[[[[[{{[[[[{\8e\08{[[[[[[{{}[[[[[{{nm,A->nm,m,M{{[{[[[{"
wait
mkdir eventfiles
q_arr=$(python $current_path/tq_gen.py -e $1)
for i in $q_arr
do
	ulim=5.0
	cd results/
	rm -rf 2_2/
	sed -i "s/#Session_number.*/#Session_number 1/g" session.dat
	echo "Setting photon energy: $i"
	source ../bin/set_momenta $1 $i
	echo "Running Initial Vegas Calculation..."
	xsval=$(source ../bin/run_vegas 5 100000 5 100000 &)
	xs=$(echo $(cut -d' ' -f1 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
	un=$(echo $(cut -d' ' -f2 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
	count=1
	if [[ $(echo $un'>'$lim1 | bc -l) == "1" ]]; then
		echo "Warning! Error not within limit, trying again..."
	fi
	while [[ $(echo $un'>'$lim1 | bc -l) == "1" && $(echo $un'>'$lim2 | bc -l) == "1" ]];
	do
		xsval=$(source ../bin/run_vegas 5 100000 5 100000 &)
		xs=$(echo $(cut -d' ' -f1 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
		un=$(echo $(cut -d' ' -f2 <<<$xsval) | sed 's/E/*10^/g;s/+//g')
		if [[ $(echo $un'<'$ulim | bc -l) == "1" ]]; then
			ulim=$un
		fi
                if [[ "$count" == "6" ]]; then
			lim2=$(echo $ulim'+'1.0 | bc -l)
                        echo "Exceeded max #iterations: $count! Trying best available upper limit: $un, $ulim , $lim2"
	        fi
                count=$(echo $count'+'1 | bc -l)
	done
        echo "Finished within error tolerance, final error: $un"
	echo "Starting event generation..."
        source ../bin/subproc_cycle $nevents > output.dat &
        wait
        mv output.dat ../eventfiles/xs_$i.txt
        mv *.lhe ../eventfiles/evt_$i.lhe
        echo "Moving results to directory complete!"
        echo "Finished calculation for neutrino energy: $1 "
	cd ..
done	
echo "Done operation!"

