#!/bin/python3

'''
Author: Sourav Sarkar
Date: October 6, 2020
Objective: Call the ranged weight factor function to store nuclear form factor
	term values in a textfile for late use in the interpolation function
'''

import numpy as np
import h5py
import EPA_XS as epa
import scipy.stats as st
from glob import glob
from optparse import OptionParser
from numpy.random import default_rng
import scipy.interpolate as ip


#define the range of photon energies
qmin=0.0449/4./1e8#GeV
#qmax=1.#GeV

#Update: increased photon energy range due to lower limit increase in the photon pdf for DIS
qmax=1.3#GeV
#number of points
print ('Generating form factor weights for Hydrogen:')
npoints=10000
#nucleus type
ntype='H'

q2,val=epa.ranged_WF(qmin,qmax,ntype,npoints)

q2=q2.reshape(-1,1)
val=val.reshape(-1,1,)

data=np.hstack((q2,val))
np.savetxt(ntype+'_weight_factor_v1.txt', data)
#nucleus type
print ('Generating form factor weights for Oxygen:')
ntype='OX'

q2,val=epa.ranged_WF(qmin,qmax,ntype,npoints)

q2=q2.reshape(-1,1)
val=val.reshape(-1,1,)

data=np.hstack((q2,val))
np.savetxt(ntype+'_weight_factor_v1.txt', data)
