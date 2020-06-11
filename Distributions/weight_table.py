#!/bin/python

'''
Author: Sourav Sarkar
'''

import numpy as np
from glob import glob
import scipy.interpolate as ip
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)
import matplotlib.colors as colors

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


datafiles=glob("Enu_*.txt")

sep_arr=np.array([])
len_arr=np.array([])
enu_arr=np.array([])

prob_arr=np.array([])
int_arr=np.array([])


def extract_array(x,y):
	rep_idx=np.where(x==x[0])
	yy=y[rep_idx]
	xx=x[:int(len(x)/len(yy))]
	dx=np.average(np.diff(xx))
	dy=np.average(np.diff(yy))
	return (dx,dy)

for i in datafiles:
	data=np.loadtxt(i,delimiter=',')
	dx,dy=extract_array(data[:,0],data[:,1])
	norm=dx*dy*np.sum(data[:,2])
	prob_norm=data[:,2]/norm

	filt=np.where(prob_norm>0)
	filt_sep=data[:,0][filt]
	filt_len=data[:,1][filt]
	filt_prob=prob_norm[filt]
	filt_enu=np.repeat(np.float_(i[4:-4]),len(filt_prob))

	sep_arr=np.append(sep_arr,filt_sep)
	len_arr=np.append(len_arr,filt_len)
	enu_arr=np.append(enu_arr,filt_enu)
	prob_arr=np.append(prob_arr,filt_prob)


sepc=sep_arr[:,np.newaxis]
lenc=len_arr[:,np.newaxis]
enuc=enu_arr[:,np.newaxis]
lenuc=np.log10(enuc)

cord=np.hstack((sepc,lenc,lenuc))
logp=np.log10(prob_arr)



############### Neutrino energy specific sampling

def weight_calculator(nu_en):
	ds=5.0
	dl=0.1
	lmax=range_func(nu_en)
	smax=np.pi*lmax*1e3
	interp_s=np.arange(0.,smax,ds)
	interp_l=np.arange(0.,lmax,dl)
	if len(interp_s)<50:
		interp_s=np.linspace(0.,smax,50)
	if len(interp_l)<50:
		interp_l=np.linspace(0.,lmax,50)

#dims=len(interp_s)
#diml=len(interp_l)
	sdiff=np.average(np.diff(interp_s))
	ldiff=np.average(np.diff(interp_l))

	smesh,lmesh=np.meshgrid(interp_s,interp_l)
	fs=smesh.flatten()[:,np.newaxis]
	fl=lmesh.flatten()[:,np.newaxis]
	fe=np.repeat(np.log10(nu_en),len(fs))[:,np.newaxis]
	interp_points=np.hstack((fs,fl,fe))

#Interpolation
	interp_lval=ip.griddata(cord,logp,interp_points,method='linear',fill_value=-np.inf)
	interp_val=10**interp_lval

###############phase space cut
	sepcut=10.#19.
	min_ang=0.01#0.019

	tot_prob=sdiff*ldiff*np.sum(interp_val)

############### detectable weight calculation
#filter 1
	filt1_id=np.where(fs<sepcut)
#filter 2
	len_cut=fs/min_ang*1e-3
	len_diff=len_cut-fl
	filt2_id=np.where(len_diff<0.)
#Adding both filters
	full_filt_id=np.union1d(filt1_id,filt2_id)
#applying the filters
	filt_interp_val=interp_val
	filt_interp_val[full_filt_id]=0.0

	det_prob=sdiff*ldiff*np.sum(filt_interp_val)

	weight=det_prob/tot_prob
#print ("Weight factor: ", tot_prob,det_prob,weight)
	return weight

en_table=np.logspace(2,7.99,250)
wt_table=np.array([])
for i in en_table:
	print ("Working on neutrino energy: ", i)
	w=weight_calculator(i)
	wt_table=np.append(wt_table,w)
enc=en_table[:,np.newaxis]
wtc=wt_table[:,np.newaxis]
data=np.hstack((enc,wtc))

np.savetxt("weight_detetion_10m.txt",data,delimiter=",")
