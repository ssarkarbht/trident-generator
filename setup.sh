#!/bin/sh

#Run this script to set up CalcHEP generator and check for all the required python modules present in the system and then run the test scripts for checking if the python module works


#Setting some print colors
green='\033[0;32m'
red='\033[0;31m'

#set basic system paths
export NTP_HOME=$(pwd)
export CALCHEP_HOME=$(pwd)/CalcHEP_generator
export CALCHEP_BUILD=$CALCHEP_HOME/calchep_3.7.1
export PYTHON_MODULE=$(pwd)/python_modules

#Set up the CalcHEP generator
echo "${green}Started setting up CalcHEP environment..."

#Set the name of the CalcHEP workspace directory
chep_workspace=NTP_CHDF
export CALCHEP_WORK=$CALCHEP_HOME/$chep_workspace

#build CalcHEP
cd $CALCHEP_HOME
tar -xvzf calchep_3.7.1.tgz
cd $CALCHEP_BUILD
make
wait
./mkWORKdir ../$chep_workspace

#Copy external SM model files
cp $CALCHEP_HOME/model_ex/* models/

echo "${green}Finished setting up CalcHEP environment..."

#Set the Neutrino trident working model and interaction process
echo "${green}Importing external Trident SM model to CalcHEP..."
cd $CALCHEP_WORK
#./calchep -blind "[[[[[[[[[[[[[[{{[[[[{\8e\08{[[[[[[{{}[[[[[{{nm,A->nm,m,M{{[{[[[{"
#./calchep -blind "[[[[[{{{\8e\08{[[[[[[{{}[[[[[{{nm,A->nm,m,M{{[{[{[}[[[{"
wait
echo "${green}Finshed setting up NTP EPA process framework..."
echo "${green}Run CalcHEP scripts for Cross-section calculation and event generation"

