#!/bin/python

#Author: Sourav Sarkar
#Date: Feb 7, 2020
# This code reads the nutrino cross section data
# from neutrino generator project in Icecube Simulation
# metaproject and make the plot of interpolated functions
# for a given energy range

import numpy as np
import matplotlib.pyplot as plt
import glob as glob
import argparse
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser("Neutrino energy")
parser.add_argument('--max_energy', '-logemax', type=float, help="Log10 of maximum of the energy range")
parser.add_argument('--min_energy', '-logemin', type=float, help="Log10 of minimum of the energy range")

args=parser.parse_args()

CC_xs_data=np.loadtxt('/home/sourav/I3/Simulation/trunk/src/neutrino-generator\
/resources/cross_section_data/csms/total_nu_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)
NC_xs_data=np.loadtxt('/home/sourav/I3/Simulation/trunk/src/neutrino-generator\
/resources/cross_section_data/csms/total_nu_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)
CCbar_xs_data = np.loadtxt('/home/sourav/I3/Simulation/trunk/src/neutrino-generator\
/resources/cross_section_data/csms/total_nubar_CC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)
NCbar_xs_data = np.loadtxt('/home/sourav/I3/Simulation/trunk/src/neutrino-generator\
/resources/cross_section_data/csms/total_nubar_NC_iso_NLO_HERAPDF1.5NLO_EIG.dat',skiprows=1)

#tot_xs_data = CC_xs_data+NC_xs_data
#totbar_xs_data = CCbar_xs_data+NCbar_xs_data



