#!/bin/bash

#clear existing environment
#module --force purge

echo "1. Loading cvmfs environment..."
#Load cvmfs environment
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

#get the job number
job_num=$1

#set path for directories
script_dir=/home/ssarkar/NTPGenerator/Job_config
mg_dir=/home/ssarkar/NTPGenerator/MadGraph/MG5_aMC_v2_8_2
p8_dir=/home/ssarkar/NTPGenerator/PYTHIA/pythia8303/examples
dataio_dir=/data/icecube/ssarkar/dataio/numu_dataset_01/batch_1
batch_dir=/data/icecube/ssarkar/madgraph_jobs/batch_1
job_dir=job_$job_num
config_file=$script_dir/numu_index1_out.h5

echo "2. Preparing run directory..."
mkdir $batch_dir/$job_dir
cd $batch_dir/$job_dir
cp -r $mg_dir/numu_* .
#cd $mg_dir

#get the parameters for madgraph run
echo "3. Reading run parameters..."
mapfile -t ev_arr < <(python $script_dir/madgraph_evt.py -f $config_file -i $job_num -n 100)
mapfile -t en_arr < <(python $script_dir/madgraph_nuen.py -f $config_file -i $job_num -n 100)
mapfile -t nc_arr < <(python $script_dir/madgraph_nucleus.py -f $config_file -i $job_num -n 100)
mapfile -t pq_arr < <(python $script_dir/madgraph_qkph.py -f $config_file -i $job_num -n 100)

#prun=1
#qrun=1
echo "4. Starting event generation..."
for i in ${!ev_arr[@]}
do
	ev=${ev_arr[$i]}
	en=${en_arr[$i]}
	nc=${nc_arr[$i]}
	pq=${pq_arr[$i]}
	echo $ev $en $pq $nc
	cd $pq
	sed -i "s/set ebeam1.*/set ebeam1 $en/g" run_script.txt
	sed -i "s/set lhaid.*/set lhaid $nc/g" run_script.txt
#	cat run_script.txt
#	if [[ $pq == 'numu_photon' ]]; then
#		prun=$(($prun + 1))
#	fi
#	if [[ $pq == 'numu_quark' ]]; then
#		qrun=$(($qrun + 1))
#	fi
	./bin/madevent run_script.txt
	wait
	cp Events/run_01/run_01_tag_1_banner.txt $dataio_dir/p8_input/event_id_$ev.txt
	cp Events/run_01/unweighted_events.lhe.gz $dataio_dir/p8_input/event_id_$ev.lhe.gz
	rm -rf Events/run_01/
	cd $dataio_dir/p8_input/
	gzip -d event_id_$ev.lhe.gz
	cd $p8_dir
	./main100 $dataio_dir/p8_input/event_id_$ev.lhe $dataio_dir/event_id_$ev.hepmc
	wait
	rm $dataio_dir/p8_input/event_id_$ev.lhe
#	echo $prun $qrun
	cd $batch_dir/$job_dir
#	cd $mg_dir
done
rm -rf $batch_dir/$job_dir/
echo "Job Successfully completed!"
