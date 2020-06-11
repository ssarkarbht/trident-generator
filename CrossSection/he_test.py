#!/bin/python3
'''
Author: Sourav Sarkar
Date: 26th Feb, 2020
Objective:
This script performs the neutrino trident coherent+diffractive cross-section
calculation at low energy where we can compare the results with some of the 
reported cross-section in the context of neutrino beam experiements
'''
#Import custom modules from trident modules directory
import sys
sys.path.insert(1,'/home/sourav/Tridents/trident_modules')
import EPA_XS as tm

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)
from glob import glob
import scipy.interpolate as ip
from scipy.integrate import quad

#Neutrino energy array
inen=np.logspace(1,8,100)

xsval=np.array([])
erval=np.array([])
unval=np.array([])
smin=0.044944
for e in inen:
	smax=2.*e*1.
	x,xerr=quad(tm.FF_Int,smin,smax,args=(e,'OX',))
	xsval=np.append(xsval,x)
	erval=np.append(erval,xerr)

#	ux,uerr=quad(tm.FF_UERR_Int,smin,smax,args=(e,'OX',))
#	lx,lerr=quad(tm.FF_LERR_Int,smin,smax,args=(e,'OX',))
#	ferr=abs(ux-lx)
#	print (e,xerr/x*100.,ferr/x*100.)
#	unval=np.append(unval,ferr)
box_data=np.loadtxt('/home/sourav/Tridents/paper_xs_results/beacom_chdf.csv', delimiter=',')
fox_data=np.loadtxt('/home/sourav/Tridents/paper_xs_results/shaofeng_ox.csv', delimiter=',')
modxs=xsval/inen*1e-36
box_func=ip.interp1d(np.log10(inen),np.log10(modxs),kind='cubic')

new_cm=np.logspace(np.log10(min(inen)),np.log10(max(inen)),1000)


f, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(7,7), gridspec_kw={'height_ratios': [4, 1]})
ax1.scatter(inen,modxs,s=5,label='Sourav')
#plt.errorbar(inen,xsval/inen*1e-36,erval/inen*1e-36,fmt='none')
#plt.errorbar(inen,xsval/inen*1e-36,unval/inen*1e-36,fmt='none')
ax1.plot(new_cm,10**(box_func(np.log10(new_cm))))
ax1.scatter(box_data[:,0],box_data[:,1],s=5,label='Zhou, Beacom')
ax1.plot(fox_data[:,0],fox_data[:,1]/fox_data[:,0]*1e-36,'r',label='Shao-Feng')
ax1.legend()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$\sigma_{\nu A}/E_{\nu}$ [$cm^2GeV^{-1}$]')
ax1.set_xlabel('Neutrino Energy (GeV)')
ax1.set_ylim([1e-44,1e-39])
ax1.grid(which='both',alpha=0.3,ls=':')

ratio=box_data[:,1]/10**(box_func(np.log10(box_data[:,0])))
ax2.scatter(box_data[:,0],ratio,s=5,c='k')
ax2.set_yscale('log')
ax2.set_ylim([0.5*min(ratio),1.5*max(ratio)])
ax2.set_ylabel('Ratio (Beacom/Sourav)')
ax2.grid(which='both',alpha=0.3,ls="--")


plt.savefig('HE_compare_ox.pdf')
