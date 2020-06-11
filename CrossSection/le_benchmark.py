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
inen=np.linspace(0.5,150,100)

xsval=np.array([])
#erval=np.array([])
cxsval=np.array([])
dxsval=np.array([])
smin=0.044944
for e in inen:
#	Z=18.
#	A=40.
#	Q2lim=(0.217/A**(1/3.))**2
	#slim=2.*e*np.sqrt(Q2lim)
	smax=2.*e*1.#(0.217/40.**(1/3.))#1.

	x,xerr=quad(tm.FF_Int,smin,smax,args=(e,'LAr',))
	xsval=np.append(xsval,x)

#	cx,cerr=quad(tm.FFC,smin,smax,args=(e,'LAr',))
#	dx,derr=quad(tm.FFD,slim,smax,args=(e,'LAr',))
        #For FI
#	ch = Z**2/137./np.pi
        #For Fdip
#	df = Z/137./np.pi
#	cxsval=np.append(cxsval,ch*cx)
#	dxsval=np.append(dxsval,df*dx)
#	val=ch*cx+df*dx
#	print (cx,dx)
#	print (x)
#	erval=np.append(erval,xerr)

#Load the matheus ntp LAr xs result
mlar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/matheus_lar_ntpxs.csv", delimiter=',')
plar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/plastid_lar_ntpxs.csv", delimiter=',')
#plar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/plestid_ar.csv", delimiter=',')
#plar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/plestid_pb.csv", delimiter=',')
slar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/Sanch_ar1.csv", delimiter=',')
#slar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/Sanch_ox.csv", delimiter=',')

#### Analytical Cross-section ######
G=1.66e-5
alpha=1/137.
Z=82.
CA=0.5
CV=0.5+2*0.2223
smax=inen*2.*1.
smin=0.0449
axs=(CA**2+CV**2)*Z**2*alpha**2*G**2/9./np.pi**3*smax*np.log10(smax/smin)


#print (erval)
#print (min(mlar_data[:,1])*1e36,max(mlar_data[:,1])*1e36)
#print (min(xs1val), max(xs1val))
plt.figure()
Z=18.#for argon
#plt.plot(inen,xsval/inen/Z**2*1e9,lw=1,linestyle='--',label="Sourav")
#plt.scatter(inen,cxsval,s=5,label="CHXS")
#plt.scatter(inen,dxsval,s=5,label="DFXS")
#plt.errorbar(inen,xsval,erval,fmt='none')
plt.plot(mlar_data[:,0],mlar_data[:,1]*1e36,lw=1,label='Matheus')
plt.plot(plar_data[:,0],plar_data[:,1]*plar_data[:,0]*Z**2*1e-9,lw=1,label='Plastid(PRD.95.073004)')
plt.plot(slar_data[:,0],slar_data[:,1]*1e36*slar_data[:,0],lw=1,label='Sanchez')
plt.plot(inen,xsval,lw=1,linestyle='--',label="Sourav")
#plt.plot(inen,axs*0.389379e9/inen/Z**2*1e9,lw=1,linestyle='-.',label="Analytical Approx.")

plt.legend()
plt.grid(which='both',alpha=0.3, linestyle=':')
plt.yscale('log')
#plt.xscale('log')
#plt.ylim([1e-4,1e2])
plt.xlim([0,20.])
plt.xlabel("Neutrino Energy (GeV)")
#plt.ylabel(r'$\sigma/(E_{\nu}Z^2)\quad [10^{-45}cm^2/GeV]$')
plt.ylabel(r"$\sigma$ [pb]")
plt.title(r"$\nu_{\mu}+LAr\rightarrow \nu_{\mu}+\mu^++\mu^-+LAr$ (Coherent+Diffractive) Cross section")
#plt.title(r"$\nu_{\mu}+Pb\rightarrow \nu_{\mu}+\mu^++\mu^-+Pb$ (Coherent+Diffractive) Cross section")
plt.savefig('LE_CHDF.pdf')
#plt.savefig('Ar_comparision1.pdf')
