#!/bin/python
'''
Author:Sourav Sarkar
Date: March 2, 2020
Objective: Check the overall weight factor that comes before the 
	   integration over CM energy (incluiding the constant
	   factor outside the integral). And plot the weight as
	   function of the photon momentum transfer for specific
	   neutrino energy.

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)
import scipy.interpolate as ip
from glob import glob
from scipy.integrate import quad
#Import custom modules from trident modules directory
import sys
sys.path.insert(1,'/home/sourav/Tridents/trident_modules')
import EPA_XS as tm


#Nu_En=10.#GeV
en=np.logspace(0,2,5)
#en=[10.]
plt.figure()
for Nu_En in en:
	print (Nu_En)
	q,val=tm.WF(Nu_En,'LAr')
	plt.plot(q,val,label='Quad %.2e' %Nu_En)
	q1,val1=tm.vegas_WF(Nu_En,'LAr')
	plt.plot(q1,val1,label='Vegas %.2e' %Nu_En)

plt.yscale('log')
plt.xscale('log')
plt.legend()
#plt.ylim([1e-2,1e2])
#plt.xlim([0.0,0.04])
plt.savefig('wf_check.pdf')
