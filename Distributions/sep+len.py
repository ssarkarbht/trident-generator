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

pltloc='/home/ssarkar/public_html/trident_plots/sep+len_plots/'
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


###################### Muon track length in Ice#################
muon_info = np.genfromtxt('mue_water_ice.txt', skip_header=10, usecols=(0,2,7,8))
E_real = muon_info[:,0] # MeV
ranges = muon_info[:,3] # cm
#Delete the elements with NaN values
E_real = np.delete(E_real,[39,98])
ranges=np.delete(ranges,[39,98])
#Convert into convenient units
E_real = E_real*1e-3 #in GeV
ranges = ranges*1e-5 #in Km
#Make the interpolation function to sample muon track length from
he_range = ip.interp1d(E_real,ranges,kind='nearest',bounds_error=False,fill_value='extrapolate')
dummy_e=np.linspace(1.01e6,1e8,10000)
dummy_r=he_range(dummy_e)
E_real=np.append(E_real,dummy_e)
ranges=np.append(ranges,dummy_r)
range_func=ip.interp1d(E_real,ranges,kind='cubic')
#################################################################


df=h5py.File("../event_extract/En_"+enustr,"r")

#Set bin edges
nsep=51
nlen=51

#set the limit bounds for the bin contruction
len_min=0.0
len_max=range_func(enu)
sep_min=0.0
sep_max=np.pi*len_max*1e3

#using np.arange instead to keep the bin size constant
#over any neutrino energies
dlen=0.1#km
dsep=5.0

lbins=np.arange(len_min,len_max,dlen)
sbins=np.arange(sep_min,sep_max,dsep)
#abins=np.linspace(0.,np.pi,nangle)#5 degree bins
#ebins=np.logspace(-1,np.log10(enu),nemin)
if len(lbins)<50:
	lbins=np.linspace(len_min,len_max,nlen)
if len(sbins)<50:
	sbins=np.linspace(sep_min,sep_max,nsep)
#Set bin centers
sbincen=0.5*(sbins[:-1]+sbins[1:])
lbincen=0.5*(lbins[:-1]+lbins[1:])

#Make the meshgrid and flatten the coordinates
smesh,lmesh=np.meshgrid(sbincen,lbincen)
sgrid=smesh.flatten()
lgrid=lmesh.flatten()

#Initialize blank arrays to store histogrammed values
q_cor=np.array([])
s_cor=np.array([])
l_cor=np.array([])
v_cor=np.array([])
#q_arr=np.array([])
for qval in df:
	#load the data
	mm=df[qval][:,0]
	mp=df[qval][:,1]
	op=df[qval][:,2]
	me=np.minimum(mm,mp)

	tracklen=range_func(me)
	seplen=op*tracklen*1e3 #convert it into meter
	#2d histogramming
	h,sepedge,lenedge=np.histogram2d(seplen,tracklen,bins=[sbins,lbins],normed=True)#,density=True)
	#flip the h matrix to have rows for range of x values with fixed y values
	h=h.T
	#Now flatten the matrix
	hflat=h.flatten()

	#filter out the zero bins (for the purpose of logscale and interpolation)
	new_h=hflat[np.where(hflat>0)]
	new_s=sgrid[np.where(hflat>0)]
	new_l=lgrid[np.where(hflat>0)]
	new_q=np.repeat(np.float_(qval),len(new_h))

	q_cor=np.append(q_cor,new_q)
	s_cor=np.append(s_cor,new_s)
	l_cor=np.append(l_cor,new_l)
	v_cor=np.append(v_cor,new_h)

#Setting up the points for grid interpolation
qc=q_cor[:,np.newaxis]
sc=s_cor[:,np.newaxis]
lc=l_cor[:,np.newaxis]

#convert q to log values for better interpolation
lqc=np.log10(qc)

cord=np.hstack((sc,lc,lqc))
log_h=np.log10(v_cor)

#setting the interpolation sampling points
#number of q grid points
nqgrid=100
#number of theta grid points
dsgrid=5.0
#number of minimum muon energy points
dlgrid=0.1

#center values where interpolation calculation is performed
interp_lq=np.linspace(min(lqc),max(lqc),nqgrid)
interp_s=np.arange(min(sc),max(sc),dsgrid)
interp_l=np.arange(min(lc),max(lc),dlgrid)

nsgrid=len(interp_s)
nlgrid=len(interp_l)
if nlgrid<50:
	nlgrid=50
	interp_l=np.linspace(min(lc),max(lc),nlgrid)
if nsgrid<50:
	nsgrid=50
	interp_s=np.linspace(min(sc),max(sc),nsgrid)
#evaluate q boundaries to calculate dq
lqedge=0.5*(interp_lq[:-1]+interp_lq[1:])
lqdiff=np.average(np.diff(interp_lq))
lqedge=np.append(lqedge,interp_lq[-1]+0.5*lqdiff)
lqedge=np.append(interp_lq[0]-0.5*lqdiff,lqedge)
interp_qdiff=np.diff(10**lqedge)

#evalute angle axis edges for 2d distribution plotting
sedge=0.5*(interp_s[:-1]+interp_s[1:])
sdiff=np.average(np.diff(interp_s))
sedge=np.append(sedge,interp_s[-1]+0.5*sdiff)
sedge=np.append(interp_s[0]-0.5*sdiff,sedge)

#evaluate min. energy axis edges for 2d dist. plotting
ledge=0.5*(interp_l[:-1]+interp_l[1:])
ldiff=np.average(np.diff(interp_l))
ledge=np.append(ledge,interp_l[-1]+0.5*ldiff)
ledge=np.append(interp_l[0]-0.5*ldiff,ledge)


ly,lq,lx=np.meshgrid(interp_l,interp_lq,interp_s)
lxc=lx.flatten()[:,np.newaxis]
lyc=ly.flatten()[:,np.newaxis]
lqc=lq.flatten()[:,np.newaxis]
interp_points=np.hstack((lxc,lyc,lqc))

#Interpolation
interp_lval=ip.griddata(cord,log_h,interp_points,method='linear',fill_value=-np.inf)
interp_mesh=interp_lval.reshape(nqgrid,nsgrid*nlgrid)

#finalprob=np.zeros(len(interp_a))
finalprob=np.zeros(nsgrid*nlgrid)
for i in range(nqgrid):
	n=10**(interp_mesh[i])
	finalprob+=n*new_wf(10**interp_lq[i])*interp_qdiff[i]


sedgemesh,ledgemesh=np.meshgrid(sedge,ledge)
probmesh=finalprob.reshape(nlgrid,nsgrid)

plt.figure()
plt.pcolormesh(sedgemesh,ledgemesh,probmesh,norm=colors.LogNorm())
plt.xlabel("Highest Track Separation (m)")
plt.ylabel("Minimum tracklength (km)")
plt.title(r"$E_{\nu}=$ %2.0e GeV" %enu)
#plt.yscale('log')
plt.colorbar()
plt.savefig(pltloc+'sep+len_%.1e.pdf' %enu)

#flatten the arrays to store into lookup table
smesh,lmesh=np.meshgrid(interp_s,interp_l)
smesh=smesh.flatten()[:,np.newaxis]
lmesh=lmesh.flatten()[:,np.newaxis]
pmesh=probmesh.flatten()[:,np.newaxis]

data=np.hstack((smesh,lmesh,pmesh))
np.savetxt("Enu_%.1e.txt" %enu, data,delimiter=',')
