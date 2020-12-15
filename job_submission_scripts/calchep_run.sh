#!/bin/bash

#Author: Sourav Sarkar
#Date: September 12, 2020
#Objective: Automate the batch script runs with varying neutrino/photon energies

echo "1-Loading CVMFS..."
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`

batch_dir=/data/icecube/ssarkar/calchep_jobs/batch_1
source_dir=/home/ssarkar/NTPGenerator/CalcHEP/workspace
job_dir=job_$1
script_dir=/home/ssarkar/NTPGenerator/Job_config
dataio_dir=/data/icecube/ssarkar/dataio/numu_dataset_01/batch_1
config_file=$script_dir/numu_index1_out.h5

echo "2-Building run directory..."
#source $source_dir/mkWORKdir $batch_dir/$job_dir
mkdir $batch_dir/$job_dir
cd $batch_dir/$job_dir
cp -r $source_dir/* .


#./calchep -blind "[[[[[{{[[[[{\8e\08{[[[[[[{{}\8a"
#cp $batch_dir/batchfile .
echo "3. Reading run parameters..."

mapfile -t i_arr < <(python $script_dir/calchep_evt.py -f $config_file -i $1 -n 100)
mapfile -t q_arr < <(python $script_dir/calchep_phen.py -f $config_file -i $1 -n 100)
mapfile -t e_arr < <(python $script_dir/calchep_nuen.py -f $config_file -i $1 -n 100)
echo "4. Starting event generation..."

for i in ${!q_arr[@]}
do
	ev=${i_arr[$i]}
        e=${e_arr[$i]}
        q=${q_arr[$i]}
        sed -i "s/p1:.*/p1:    $e/g" batchfile
        sed -i "s/p2:.*/p2:    $q/g" batchfile
        echo "Running energies $e $q  ..."
        ./calchep_batch batchfile
        wait
        echo "Done calculation! Moving the results to database..."
        cp batch_results/*.lhe.gz $dataio_dir/event_id_$ev.lhe.gz
	cp html/runs/single.txt $dataio_dir/p8_input/event_id_$ev.txt
	
done
rm -rf $batch_dir/$job_dir/
echo "Done operation!"

