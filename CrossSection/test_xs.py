#!/bin/python

#Author: Sourav Sarkar
# Feb 21, 2020
#Desciption: The purpose of this script is to provide required modules
#that runs appropriate functions and integation method to evaluate the 
#final incoming neutrino-nucleon cross-section using the EPA model
#for coherent and diffractive regime

#User should have the option to choose different form factor models
#

import numpy as np
import scipy.interpolate as ip
from scipy.integrate import dblquad
from glob import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)


#################################---Segment for building EPA cross-section function----#########
#Loading the CM EPA cross-section data
epa_data=np.loadtxt("/home/sourav/Tridents/new_trident_mdl/epa_data.dat")

cmen=epa_data[:,0]
cmxs=epa_data[:,1]

#Defining the interpolation function in log scale (both in energy and cross-section)
epa_func=ip.UnivariateSpline(np.log10(cmen),np.log10(cmxs),s=0)

#define a wrapper function to spit out cross-section in pb
def epa_xs(s):
        lcmen=np.log10(np.sqrt(s))
        val=epa_func(lcmen)
        return 10**val

#################################################################################################
print (epa_xs(1.0))
'''
dumm = np.logspace(0,3,1000)
plt.figure()
plt.scatter(dumm,epa_xs(dumm**2),s=3)
plt.savefig('check.pdf')
'''

##########################-----Form Factors------########################################
#Making the nuclear type library for form factor calculation
ntype_lib={'OX':(16.,8.),'H':(1.,1.), 'LAr':(40.,18.)}


#Coherent Form Factor:
#// Phys.Rev. D3 (1971) 2686-2706

#The function takes the photon momentum transfer squared and nucleus type
#Implemented Nucleus type is Oxygen and Hydrogen for now (for the purpose of IceCube Detector)
def FI_I(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	a = (5.076*(0.58+0.82*A**(1/3.)))**2
	fac = (1+a*Q2/12.)**2
	return (1/fac)**2/Q2

#Wood-Saxon form factor in the analytical form from Matheus Hoestert article
#
def FI_II(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	a = 0.523*5.068 #1fm=5.068GeV^-1
	r0 = 1.126*A**(1/3.)*5.068
	q=np.sqrt(Q2)
	val=(3*np.pi*a/(r0**2+np.pi*np.pi*a**2))*(np.pi*a*(1/np.tanh(np.pi*q*a))*np.sin(q*r0)-r0*np.cos(q*r0))/(q*r0*np.sinh(np.pi*q*a))
	return val**2/Q2


'''
q=np.logspace(-9,1,1000)
plt.figure()
plt.plot(q,q**2*FI_II(q**2,'OX'))
plt.plot(q,q**2*FI_I(q**2,'OX'))
#plt.plot(q,FI_I(q**2,'OX')/FI_I(q[0]**2,'OX'))
plt.yscale('log')
plt.xscale('log')
#plt.xlim([1e-2,7e-1])
plt.ylim([1e-5,1e2])
plt.savefig('testing.pdf')
'''

#Diffractive dipole form factor
#// 1612.05642

def Fdip(Q2,ntype):
	ninfo=ntype_lib.get(ntype)
	if not(ninfo):
		print ('Error: No know nucleus type ',ntype,' found in the library')
		print ('Available Nucleus types are: ', ntype_lib.keys())
		return 0
	else:
		A,Z=ninfo[0], ninfo[1]
	gdip = (1+Q2/0.71)**(-2)
	tau  = Q2/4.
	xi   = 4.7
	fac  = gdip*(1+tau*xi)/(1+tau)
	return fac**2/Q2

##############################################################################################

def Tot_ch(Q2,s,en,ntype):
	return epa_xs(s)*FI_II(Q2,ntype)/s
def Tot_df(Q2,s,en,ntype):
	return epa_xs(s)*Fdip(Q2,ntype)/s

def Q2chl(s):
	global e
	return (s/2./e)**2
def Q2chu(s):
	global A
	return (0.217/A**(1/3.))**2
def Q2dfl(s):
	global e, A
	temp1=(0.217/A**(1/3.))**2
	temp2=(s/2./e)**2
	return max(temp1,temp2)
def Q2dfu(s):
	return 1.69
#Neutrino energy array
inen=np.linspace(0.5,15,30)

xsval=np.array([])
#erval=np.array([])
cxsval=np.array([])
dxsval=np.array([])
smin=0.044944
for e in inen:
        Z=18.
        A=40.
        Q2lim=(0.217/A**(1/3.))**2
        slim=2.*e*np.sqrt(Q2lim)
        smax=2.*e*1.3#(0.217/40.**(1/3.))#1.

        cx,cerr=dblquad(Tot_ch,smin,slim,Q2chl,Q2chu,args=(e,'LAr',))
        dx,derr=dblquad(Tot_df,smin,smax,Q2dfl,Q2dfu,args=(e,'LAr',))
        #For FI
        ch = Z**2/137./np.pi
        #For Fdip
        df = Z/137./np.pi
        cxsval=np.append(cxsval,ch*cx)
        dxsval=np.append(dxsval,df*dx)

#Load the matheus ntp LAr xs result
mlar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/matheus_lar_ntpxs.csv", delimiter=',')
plar_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/plastid_lar_ntpxs.csv", delimiter=',')

#print (erval)
#print (min(mlar_data[:,1])*1e36,max(mlar_data[:,1])*1e36)
#print (min(xs1val), max(xs1val))
plt.figure()
#plt.scatter(inen,xsval,s=5,label="Sourav1")
plt.scatter(inen,cxsval,s=5,label="CHXS")
plt.scatter(inen,dxsval,s=5,label="DFXS")
#plt.errorbar(inen,xsval,erval,fmt='none')
plt.scatter(mlar_data[:,0],mlar_data[:,1]*1e36,s=5,label='Matheus')
plt.scatter(plar_data[:,0],plar_data[:,1]*1e-9*plar_data[:,0]*18.**2,s=5,label='Plastid')

plt.legend()
plt.grid(which='both',alpha=0.3)
plt.yscale('log')
plt.ylim([1e-9,1e-1])
plt.xlim([0,15.5])
plt.savefig('test.pdf')

