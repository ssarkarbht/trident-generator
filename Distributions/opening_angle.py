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

plt.figure()
plt.scatter(cmen**2/4./enu,epa_xs(cmen**2),s=5)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-15,1e-2])
plt.xlim([0.0449/4./enu,1.])
plt.savefig(pltloc+"xs_check.pdf")


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
#df=h5py.File("../event_extract/En_10.00","r")

#Set bin edges
abins=np.linspace(0.,np.pi,37)#5 degree bins
#ebins=np.logspace(-1,np.log(enu),51)

abincen=0.5*(abins[:-1]+abins[1:])
#ebincen=0.5*(ebins[:-1]+ebins[1:])

q_cor=np.array([])
a_cor=np.array([])
v_cor=np.array([])
q_arr=np.array([])
print (len(df))
#def colors(q):
#	global enu, qmin
#	qdiff=0.0-np.log10(qmin)
#	val= (np.log10(np.float_(q))-np.log10(qmin))/qdiff
#	return plt.cm.jet(val)
colors=plt.cm.jet(np.linspace(0,1,len(df)))
plt.figure()
c=0
for qval in df:
#	mm=df[qval][:,0]
#	mp=df[qval][:,1]
	op=df[qval][:,2]
#	me=np.minimum(mm,mp)

#	h,opedge,meedge=np.histogram2d(op,me,bins=[abins,ebins],density=True)
#	nabins=abincen()
	h,bins=np.histogram(op,bins=abins,density=True)

	plt.step(abincen*180./np.pi,h,where='mid',color=colors[c],label=qval)

	new_h=h[np.where(h>0)]
	new_a=abincen[np.where(h>0)]
	new_q=np.repeat(np.float_(qval),len(new_h))

	q_cor=np.append(q_cor,new_q)
	a_cor=np.append(a_cor,new_a)
	v_cor=np.append(v_cor,new_h)
	c+=1
plt.yscale('log')
plt.xlabel("Opening angle (degree)")
plt.title("Nu_E = %.1e GeV" %enu)
plt.legend(fontsize=3)
plt.savefig(pltloc+'op_check.pdf')

#Setting the points for grid interpolation
qc=q_cor[:,np.newaxis]
ac=a_cor[:,np.newaxis]
lqc=np.log10(qc)
cord=np.hstack((ac,lqc))
log_h=np.log10(v_cor)

#setting the interpolation sampling points
#number of q grid points
nqgrid=250
#number of theta grid points
nagrid=90

interp_lq=np.linspace(min(lqc),max(lqc),nqgrid)
#interp_a=np.linspace(min(ac),max(ac),nagrid)
interp_ae=np.linspace(0.,np.pi,nagrid+1)
interp_a=0.5*(interp_ae[:-1]+interp_ae[1:])

lqedge=0.5*(interp_lq[:-1]+interp_lq[1:])
lqdiff=np.average(np.diff(interp_lq))
lqedge=np.append(lqedge,interp_lq[-1]+0.5*lqdiff)
lqedge=np.append(interp_lq[0]-0.5*lqdiff,lqedge)

interp_qdiff=np.diff(10**lqedge)

lx,ly=np.meshgrid(interp_a,interp_lq)
lxf=lx.flatten()
lyf=ly.flatten()
lxc=lxf[:,np.newaxis]
lyc=lyf[:,np.newaxis]
interp_points=np.hstack((lxc,lyc))

#Interpolation
interp_lval=ip.griddata(cord,log_h,interp_points,method='linear',fill_value=-np.inf)
interp_mesh=interp_lval.reshape(nqgrid,nagrid)

finalprob=np.zeros(len(interp_a))
for i in range(len(interp_lq)):
	n=10**(interp_mesh[i])
	finalprob+=n*new_wf(10**interp_lq[i])*interp_qdiff[i]

plt.figure()
plt.step(interp_a*180./np.pi,finalprob,where='mid',c='g')
plt.yscale('log')
plt.title(r"$E_{\nu}=$ %.1e GeV" %enu)
plt.ylabel('Diff. distribution of opening angle')
plt.xlabel("Opening angle between two muon tracks (degree)")
plt.grid(which='both',alpha=0.3,ls='--')
plt.savefig(pltloc+'op_ang_%.1e.pdf' %enu)


