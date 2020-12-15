#!usr/bin/env python

import numpy as np
import h5py as h5
from optparse import OptionParser

#command line options to have input parameters
parser=OptionParser()
parser.add_option("-f", "--InputFilename", dest="INFILE", type=str)

#parse the parameters
(options,args)=parser.parse_args()

#load the hdf5 file
infile=h5.File(options.INFILE, "r")

#extract event index
data_arr=infile['DIS_interactions']['interaction_param']['event_id']

#failsafe for last job segment in case iteration ends before reaching
#total events per job value
data_len=len(data_arr)

print ("No. of DIS interaction for this batch: ", data_len)
