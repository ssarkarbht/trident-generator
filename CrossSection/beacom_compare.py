import numpy as np
import scipy.interpolate as ip
from scipy.integrate import quad
from glob import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('figure',dpi=350)


#################################---Segment for building EPA cross-section function----#########
#Loading the CM EPA cross-section data
epa_data=np.loadtxt("/home/sourav/Tridents/new_trident_mdl/epa_data_new1.dat")

cmen=epa_data[:,0]
cmxs=epa_data[:,1]

#Defining the interpolation function in log scale (both in energy and cross-section)
epa_func=ip.UnivariateSpline(np.log10(cmen),np.log10(cmxs),s=0)

#define a wrapper function to spit out cross-section in pb
def epa_xs(s):
        lcmen=np.log10(np.sqrt(s))
        val=epa_func(lcmen)
        return 10**val

beacom_data=np.loadtxt("/home/sourav/Tridents/paper_xs_results/beacom_epa.csv", delimiter=',')
dumm=np.logspace(-1,4,10000)

#res=6.5e3
res=6.46e3
#two options for different format of cross-section plotting
#op=1#if plain cross-section in y-axis
op=2#if xs/s
plt.figure()
if op==1:
	ymin=1e-47
	ymax=1e-34
	plt.plot(dumm**2,epa_xs(dumm**2)*1e-36,ls='--',c='b',label='Sourav')
	plt.scatter(beacom_data[:,0],beacom_data[:,1]*beacom_data[:,0],s=5,c='k',label='Zhou,Beacom')
	plt.vlines(res,ymin,ymax,'r',label="W-resonance")
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim([ymin,ymax])
	plt.title(r"$\nu_{\mu}+\gamma(real)\rightarrow\nu_{\mu}+\mu^++\mu^-$ EPA Comparision")
	plt.xlabel(r"Squared CM energy (s) ($GeV^2$)")
	plt.ylabel(r"$\sigma_{\nu\gamma}$ [$cm^2$]")
	plt.legend()
	plt.grid(which='both',alpha=0.3,ls=':')

elif op==2:
	ymin=1e-45
	ymax=1e-38
	f,(ax1,ax2)=plt.subplots(nrows=2, sharex=True, figsize=(7,7), gridspec_kw={'height_ratios': [4, 1]})
	ax1.plot(dumm**2,epa_xs(dumm**2)/dumm**2*1e-36,ls='--',c='b',label='Sourav')
	ax1.scatter(beacom_data[:,0],beacom_data[:,1],s=5,c='k',label='Zhou,Beacom')
	ax1.vlines(res,ymin,ymax,'r',label="W-resonance")
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax1.set_ylim([ymin,ymax])
	ax1.set_title(r"$\nu_{\mu}+\gamma(real)\rightarrow\nu_{\mu}+\mu^++\mu^-$ EPA Comparision")
	ax1.set_ylabel(r"$\sigma_{\nu\gamma}/s$ [$cm^2GeV^{-2}$]")
	ax1.set_xlabel(r"Squared CM energy (s) ($GeV^2$)")
	ax1.legend()
	ax1.grid(which='both',alpha=0.3,ls=':')

	ratio=beacom_data[:,1]/(epa_xs(beacom_data[:,0])/beacom_data[:,0]*1e-36)
	print (min(ratio),max(ratio))
	ax2.scatter(beacom_data[:,0],ratio,s=5,c='k')
	ax2.set_yscale('log')
	ax2.set_ylim([3e-2,7e0])
	ax2.set_ylabel("Ratio (Beacom/Sourav)")
	ax2.grid(which='both',alpha=0.3)
plt.savefig('epa_compare.pdf')

