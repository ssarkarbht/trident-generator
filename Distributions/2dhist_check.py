#!/bin/python3
'''
Author: Sourav Sarkar
'''

import numpy as np
from glob import glob
import scipy.interpolate as ip
import scipy.stats as st
import EPA_XS as tm
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)
import matplotlib.colors as colors


pltloc='/home/ssarkar/public_html/trident_plots/'
enu=10.00
df=h5py.File("../event_extract/En_10.00","r")

#Set bin edges
abins=np.linspace(0.,np.pi,37)#5 degree bins
ebins=np.logspace(-1,np.log10(enu),51)

abincen=0.5*(abins[:-1]+abins[1:])
ebincen=0.5*(ebins[:-1]+ebins[1:])

qval=list(df.keys())

mm=df[qval[0]][:,0]
mp=df[qval[0]][:,1]
op=df[qval[0]][:,2]
me=np.minimum(mm,mp)

print (max(me))
h,opedge,meedge=np.histogram2d(op,me,bins=[abins,ebins])#,normed=True)
h=h.T

amesh,emesh=np.meshgrid(abincen,ebincen)
amesh=amesh.flatten()[:,np.newaxis]
emesh=emesh.flatten()[:,np.newaxis]
hmesh=h.flatten()[:,np.newaxis]

data=np.hstack((amesh,emesh,hmesh))
print (data)

print (max(meedge))
x,y=np.meshgrid(opedge,meedge)
plt.figure()
plt.pcolormesh(x,y,h,norm=colors.LogNorm())
plt.yscale('log')
plt.xlabel("Opening angle (degree)")
plt.ylabel("Min. Muon Energy (GeV)")
plt.colorbar()
plt.savefig(pltloc+"2dhist_check_10.pdf")


plt.figure()
plt.scatter(abincen,np.sum(h,axis=0))
plt.hist(op,bins=abins,histtype='step')
plt.yscale('log')
plt.savefig(pltloc+'angle_check_10.pdf')

