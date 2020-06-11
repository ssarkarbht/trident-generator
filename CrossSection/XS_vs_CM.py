#!/bin/python3

'''
Author: Sourav Sarkar
Date: 24th Feb, 2020
Objective:This script loads the raw EPA cross-section data from CalcHEP
	  Then plots data from all the incoming neutrino energies 
	  makes the interpolation of the cross-section as afunction of
	  CM energy of the neutrin-photon system and stores the
	  interpolated data into a new file that can be loaded and
	  interpolated for later use down the final XS calculation steps 
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)
from glob import glob
#from scipy.interpolate import interp1d
#from scipy.interpolate import UnivariateSpline
#from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.interpolate as intp

def data_extract(files):
	q = np.array([])
	x = np.array([])
	u = np.array([])
	for i in files:
		tmpf=open(i,'r')
		tmpl=tmpf.readlines()[0].replace("*10^",'e')
		infstr = tmpl.split()
		q = np.append(q,np.float_(infstr[1]))
		x = np.append(x,np.float_(infstr[2]))
		u = np.append(u,np.float_(infstr[3]))
	return (np.float_(infstr[0]),q,x,u)

#This location has the EPA result where CaclHEP had zero muon mass (wrong approx.)
#locs = sorted(glob('/home/sourav/CalcHEP/NTP_CHDF/workspace_v0/EPA_XS1/*'))

#Updated result with corrected Muon mass in CalcHEP model
#locs = sorted(glob('/home/sourav/CalcHEP/NTP_CHDFv1/workspacev0/EPA_XS0/*'))

#Updated result with Muon mass + CM energy correction
locs = sorted(glob('/home/sourav/CalcHEP/NTP_CHDFv1/workspacev0/EPA_XS/*'))

#Store the energy and momentum trasnfer values all in single arrays
cmdict = {}
xs_arr = np.array([])
cm_arr = np.array([])
un_arr = np.array([])
#interpolation order
order='slinear'
#order='quadratic'
#order='cubic'
plt.figure()

for l in locs:
	files=glob(l+'/xs_*.txt')
	evalu=data_extract(files)
	cmdict[evalu[0]]=evalu[1]
	sort_index=np.argsort(evalu[1])
	q = evalu[1][sort_index]
	x = evalu[2][sort_index]
	u = evalu[3][sort_index]

	q=q[np.where(x!=0.)]
	x=x[np.where(x!=0.)]
	u=u[np.where(x!=0.)]

#	cm=np.sqrt(2.*q*evalu[0])
	cm=np.sqrt(4.*q*evalu[0])
	cm_arr=np.append(cm_arr,cm)
	xs_arr=np.append(xs_arr,x)
	un_arr=np.append(un_arr,u)

	plt.scatter(cm,x,s=5,c='k')#,label=str(evalu[0]))
	#plt.errorbar(cm,x,x*u/100.,c='k',fmt='none')

#Following part was used to make individual interpolated plots and compare to see if they match with each other as a function of CM energy

#	xs_func=interp1d(np.log10(evalu[1]),np.log10(evalu[2]),kind=order)
#	q1=np.logspace(np.log10(min(evalu[1])),np.log10(max(evalu[1])),1000)
#	plt.plot(np.sqrt(2.0*q1*evalu[0]),10**xs_func(np.log10(q1)),label=str(evalu[0]))
#	plt.fill_between(2.0*q*evalu[0], x*(1+u/100.), x*(1-u/100.), alpha=0.3)

#### Data Cleaning (getting rid of outliers)
def clean_data(c,x,u):
	inip=51
#	cl=np.log10(c)
#	xl=np.log10(x)
	cmavg=np.average(c[:inip])
	xsavg=np.average(x[:inip])
	unavg=1./np.sqrt(np.sum(1./u[0:inip]**2)/inip)

	t1=np.delete(c,np.arange(inip))
	t2=np.delete(x,np.arange(inip))
	t3=np.delete(u,np.arange(inip))
	cn=np.append(cmavg,t1)
	xn=np.append(xsavg,t2)
	un=np.append(unavg,t3)

	subc=cn[np.where(cn<1.0)]
	subx=xn[np.where(cn<1.0)]
	diff=(subx[:-1]-subx[1:])
	diff=np.append(diff,0.0)

	tol=0.01
	bdp=np.where(diff>tol)
	cn=np.delete(cn,bdp)
	xn=np.delete(xn,bdp)
	un=np.delete(un,bdp)
	#mean=np.append(mean,xn[-1])
	#ratio=xn/mean
#	plt.figure()
#	plt.scatter(subc,diff,s=5)
#	plt.savefig('filter_data.pdf')
	return (cn,xn,un)


#Check the resonance energy
cmres=80.62
plt.vlines(cmres,1e-9,1e2,'r',lw=1,label='W-resonance')

cmres=0.21189
#plt.vlines(cmres,1e-6,1e2,'b',lw=1)
'''
This part performs the univariate spline interpolation on the entire population of xs vs cm data
for all incoming neutrino energies
'''

sorted_index=np.argsort(cm_arr)
new_cmarr=np.log10(cm_arr[sorted_index])
new_xsarr=np.log10(xs_arr[sorted_index])
new_unarr=un_arr[sorted_index]

#print (clean_data(new_cmarr,new_xsarr))
clean=clean_data(new_cmarr,new_xsarr,new_unarr)
print (min(10**clean[0]))
plt.scatter(10**clean[0],10**clean[1],s=5,c='r')
#########################INTERPOLATION#############################
#for old incorrect dataset (Qmax=1.3GeV)
#xs_usi = intp.UnivariateSpline(clean[0],clean[1],w=clean[2],s=120)

#for new muon mass corrected dataset (Qmax=1.3GeV)
xs_usi = intp.UnivariateSpline(clean[0],clean[1],s=1)

###################################################################
cml=np.logspace(min(clean[0]),max(clean[0]),10000)
xsl=10**xs_usi(np.log10(cml))
plt.plot(cml,xsl,'b',label="Univariate Spline")


###############Analytical comparison####################
cma=np.logspace(min(clean[0]),2,25)
CA=0.5
CV=0.5+2*0.2223

G=1.166e-5
alpha=1/137.
smin=0.044944
xsa=G**2*alpha*cma**2/np.pi**2/9.*np.log(cma**2/smin)
#xsa=0.5*(CA**2+CV**2)*(2.*G**2*alpha/9./np.pi/np.pi)*cma**2*(np.log(cma**2/0.011236)-19./6)
plt.plot(cma,xsa*0.389379e9,'g',label="Analytical XS")
#Store the interpolated data for later use
col_cml = cml.reshape(len(cml),1)
col_xsl = xsl.reshape(len(cml),1)
epa_data = np.hstack((col_cml,col_xsl))
np.savetxt('epa_data_new1.dat',epa_data)

#plt.legend(fontsize=3.5,loc='upper left')
plt.xlabel(r'CM Energy $\sqrt{s}$ (GeV)')
plt.yscale('log')
plt.xscale('log')
#plt.xlim([1e-2,1e9])
plt.ylim([1e-9,1e2])
plt.title(r'Cross-section of ($\nu_{\mu}+\gamma\rightarrow\nu_{\mu}+\mu^++\mu^-$) as a function of CM energy')
#standard limit
#plt.ylim([1e-6,1e2])

#for testing purpose: limiting the range of the plot:
#plt.xlim([1e-1,1e0])
#plt.ylim([1e-6,1e-4])

#plt.xlim([1e0,1e3])
#plt.ylim([1e-4,1e-1])

plt.grid(which='both',alpha=0.3,linestyle='--')
plt.ylabel(r'EPA $\sigma$ (pb) $[10^{-36} cm^2]$')
plt.legend()
plt.savefig('CMXS_plot.pdf')
plt.close()

#Plotting the evaluated CM energy values:
nuen=np.logspace(0,8,1000)
phqm=0.0449/2./nuen
plt.figure()
for key, value in cmdict.items():
	plt.scatter(np.repeat(key,len(value)),value,s=5)
plt.plot(nuen,phqm,'r',label="CM Threshold")
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-9,2e0])
plt.xlabel('Neutrino Energy (GeV)')
plt.ylabel('Photon Momentum Transfer (GeV)')
plt.grid(which='both',alpha=0.3,ls=":")
plt.legend()
plt.savefig("Scatter_CM.pdf")



#plt.scatter(clean)
#taking average of the values for the first CM energy, as they are all the same 
#for all incoming neutrino energies
#fcm = np.average(new_cmarr[0:42])
#fxs = np.average(new_xsarr[0:42],weights=1/new_unarr[0:42])
#fun = 1./np.sqrt(np.sum(1./new_unarr[0:42]**2)/42.)
#print (fxs)

#t1=np.delete(new_cmarr,np.arange(42))
#t2=np.delete(new_xsarr,np.arange(42))
#t3=np.delete(new_unarr,np.arange(42))

#t11=np.append(fcm,t1)
#t22=np.append(fxs,t2)
#t33=np.append(fun,t3)

#respoint (in logscale)
#resp=1.778
#resp=4.0
#fil_cmarr=t11[np.where(t11<resp)]
#fil_xsarr=t22[np.where(t11<resp)]
#fil_unarr=t33[np.where(t11<resp)]

#Testing purpose:
#print (len(fil_cmarr))
#print (min(fil_cmarr),max(fil_cmarr))
#xs_usi = intp.UnivariateSpline(new_cmarr,new_xsarr,k=3)#,s=840)
#cml=np.logspace(min(new_cmarr),max(new_cmarr),1000)
#print (fil_cmarr)

