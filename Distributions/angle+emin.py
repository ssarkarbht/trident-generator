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
import argparse

parser=argparse.ArgumentParser("Neutrino Energy")
parser.add_argument('--energy', '-e', type=str)
args=parser.parse_args()

enustr = args.energy
enu=np.float_(enustr)

pltloc='/home/ssarkar/public_html/trident_plots/'
#enu=10.00
###################### Q Weight Factor ######################

#Loading the CM EPA cross-section data
epa_data=np.loadtxt("epa_data_new1.dat")

cmen=epa_data[:,0]
cmxs=epa_data[:,1]
#Defining the interpolation function in log scale (both in energy and cross-section)
epa_func=ip.UnivariateSpline(np.log10(cmen),np.log10(cmxs),s=0)
#define a wrapper function to spit out cross-section in pb
def epa_xs(s):
        lcmen=np.log10(np.sqrt(s))
        val=epa_func(lcmen)
        return 10**val

'''
plt.figure()
plt.scatter(cmen**2/4./enu,epa_xs(cmen**2),s=5)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-15,1e-2])
plt.xlim([0.0449/4./enu,1.])
plt.savefig(pltloc+"xs_check.pdf")
'''

q2,val=tm.WF(enu,'OX')
q=np.sqrt(q2)
qmin=min(q)
qmax=max(q)
weight_func=ip.interp1d(np.log10(q),np.log10(val))
def new_wf(qval):
        global enu
        scm=4.*enu*qval
        return epa_xs(scm)*10**weight_func(np.log10(qval))

'''
plt.figure()
#plt.scatter(q,epa_xs(q*4.*enu))
#plt.plot(q,10**weight_func(np.log10(q)))
plt.plot(q,new_wf(q),label='interpolation')
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-15,1e-2])
plt.xlim([0.0449/4./enu,1.])

#plt.legend()
plt.savefig(pltloc+"weight_check.pdf")
'''
###############################################################


df=h5py.File("../event_extract/En_"+enustr,"r")

#Set bin edges
nangle=37
nemin=51

abins=np.linspace(0.,np.pi,nangle)#5 degree bins
ebins=np.logspace(-1,np.log10(enu),nemin)

#Set bin centers
abincen=0.5*(abins[:-1]+abins[1:])
ebincen=0.5*(ebins[:-1]+ebins[1:])

#Make the meshgrid and flatten the coordinates
amesh,emesh=np.meshgrid(abincen,ebincen)
agrid=amesh.flatten()
egrid=emesh.flatten()

#Initialize blank arrays to store histogrammed values
q_cor=np.array([])
a_cor=np.array([])
e_cor=np.array([])
v_cor=np.array([])
#q_arr=np.array([])

for qval in df:
	#load the data
	mm=df[qval][:,0]
	mp=df[qval][:,1]
	op=df[qval][:,2]
	me=np.minimum(mm,mp)

	#2d histogramming
	h,opedge,meedge=np.histogram2d(op,me,bins=[abins,ebins],normed=True)#,density=True)
	#flip the h matrix to have rows for range of x values with fixed y values
	h=h.T
	#Now flatten the matrix
	hflat=h.flatten()

	#filter out the zero bins (for the purpose of logscale and interpolation)
	new_h=hflat[np.where(hflat>0)]
	new_a=agrid[np.where(hflat>0)]
	new_e=egrid[np.where(hflat>0)]
	new_q=np.repeat(np.float_(qval),len(new_h))

	q_cor=np.append(q_cor,new_q)
	a_cor=np.append(a_cor,new_a)
	e_cor=np.append(e_cor,new_e)
	v_cor=np.append(v_cor,new_h)

#Setting up the points for grid interpolation
qc=q_cor[:,np.newaxis]
ac=a_cor[:,np.newaxis]
ec=e_cor[:,np.newaxis]

#convert energy and q to log values for better interpolation
lqc=np.log10(qc)
lec=np.log10(ec)

cord=np.hstack((ac,lec,lqc))
log_h=np.log10(v_cor)

#setting the interpolation sampling points
#number of q grid points
nqgrid=100
#number of theta grid points
nagrid=90
#number of minimum muon energy points
negrid=100

#center values where interpolation calculation is performed
interp_lq=np.linspace(min(lqc),max(lqc),nqgrid)
interp_a=np.linspace(min(ac),max(ac),nagrid)
interp_le=np.linspace(min(lec),max(lec),negrid)

#evaluate q boundaries to calculate dq
lqedge=0.5*(interp_lq[:-1]+interp_lq[1:])
lqdiff=np.average(np.diff(interp_lq))
lqedge=np.append(lqedge,interp_lq[-1]+0.5*lqdiff)
lqedge=np.append(interp_lq[0]-0.5*lqdiff,lqedge)
interp_qdiff=np.diff(10**lqedge)

#evalute angle axis edges for 2d distribution plotting
aedge=0.5*(interp_a[:-1]+interp_a[1:])
adiff=np.average(np.diff(interp_a))
aedge=np.append(aedge,interp_a[-1]+0.5*adiff)
aedge=np.append(interp_a[0]-0.5*adiff,aedge)

#evaluate min. energy axis edges for 2d dist. plotting
leedge=0.5*(interp_le[:-1]+interp_le[1:])
lediff=np.average(np.diff(interp_le))
leedge=np.append(leedge,interp_le[-1]+0.5*lediff)
leedge=np.append(interp_le[0]-0.5*lediff,leedge)
eedge=10**leedge


ly,lq,lx=np.meshgrid(interp_le,interp_lq,interp_a)
lxc=lx.flatten()[:,np.newaxis]
lyc=ly.flatten()[:,np.newaxis]
lqc=lq.flatten()[:,np.newaxis]
interp_points=np.hstack((lxc,lyc,lqc))

#Interpolation
interp_lval=ip.griddata(cord,log_h,interp_points,method='linear',fill_value=-np.inf)
interp_mesh=interp_lval.reshape(nqgrid,nagrid*negrid)

#finalprob=np.zeros(len(interp_a))
finalprob=np.zeros(nagrid*negrid)
for i in range(nqgrid):
	n=10**(interp_mesh[i])
	finalprob+=n*new_wf(10**interp_lq[i])*interp_qdiff[i]


aedgemesh,eedgemesh=np.meshgrid(aedge,eedge)
probmesh=finalprob.reshape(negrid,nagrid)

plt.figure()
plt.pcolormesh(aedgemesh*180./np.pi,eedgemesh,probmesh,norm=colors.LogNorm())
plt.xlabel("Opening angle (Degree)")
plt.ylabel("Minimum Muon Energy (GeV)")
plt.title(r"$E_{\nu}=$ %2.0e GeV" %enu)
plt.yscale('log')
plt.colorbar()
plt.savefig(pltloc+'angle+emin_10.pdf')
