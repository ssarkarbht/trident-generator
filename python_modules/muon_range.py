#!/bin/python3
'''
Author: Sourav Sarkar
Date: 25-05-2020
Objective: Module to get muon track length as a function of muon energy
'''

import numpy as np
import scipy.interpolate as ip

#Load the PDG data"
muon_info = np.genfromtxt('mue_water_ice.txt', skip_header=10, usecols=(0,2,7,8))
#Extract Muon energies and corresponding range
E_real = muon_info[:,0] # MeV
ranges = muon_info[:,3] # cm
#Delete the elements with NaN values
E_real = np.delete(E_real,[39,98])
ranges=np.delete(ranges,[39,98])
#Convert into convenient units
E_real = E_real*1e-3 #in GeV
ranges = ranges*1e-5 #in Km
#Make the interpolation function to extrpolate muon track length at
#high energy beyond the energy range of data (should not be heavily relied on)
he_range = ip.interp1d(E_real,ranges,kind='nearest',bounds_error=False,fill_value='extrapolate')
#generate some dummy high energy hypothetical muon range data (in the energy range
# of 1e6-1e8 GeV)
dummy_e=np.linspace(1.01e6,1e8,10000)
dummy_r=he_range(dummy_e)
#add the dummy samples to rest of the data
E_real=np.append(E_real,dummy_e)
ranges=np.append(ranges,dummy_r)
#make the final cubic interpolation
range_func=ip.interp1d(E_real,ranges,kind='cubic')

