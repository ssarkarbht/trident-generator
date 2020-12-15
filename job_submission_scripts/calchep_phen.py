#!usr/bin/env python

import numpy as np
import h5py as h5
from optparse import OptionParser

#command line options to have input parameters
parser=OptionParser()
parser.add_option("-f", "--InputFilename", dest="INFILE", type=str)
parser.add_option("-i", "--StartIndex", dest="INDEX",type=int)
parser.add_option("-n", "--EventNumber", dest="NUMBER",type=int)

#parse the parameters
(options,args)=parser.parse_args()
index=options.INDEX
evtnum=options.NUMBER

#load the hdf5 file
infile=h5.File(options.INFILE, "r")

#extract photon energies
data_arr=infile['CH+DF_interactions']['interaction_param']['photon_energy']

#failsafe for last job segment in case iteration ends before reaching
#total events per job value
data_len=len(data_arr)
start_index=evtnum*index
end_index=min(start_index+evtnum-1,data_len-1)
iter_len=(end_index-start_index)+1

for i in range(iter_len):
    evt_index=start_index+i
    print (data_arr[evt_index])
