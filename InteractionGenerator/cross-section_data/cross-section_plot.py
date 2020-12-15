#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc("figure",dpi=250)
import scipy.interpolate as ip
pltloc='/home/ssarkar/public_html/trident_plots/diff_distributions/'


total_data=np.loadtxt('beacom_total.csv',delimiter=',')
disqk_data=np.loadtxt('new_beacom_DIS_quark.csv',delimiter=',')
disph_data=np.loadtxt('new_beacom_DIS_photon.csv',delimiter=',')
chdf_data=np.loadtxt('new_beacom_chdf.csv',delimiter=',')


def serialize(x,y):
	index=np.argsort(x)
	return (x[index],y[index])

totalx,totaly=serialize(total_data[:,0],total_data[:,1])
disqkx,disqky=serialize(disqk_data[:,0],disqk_data[:,1])
disphx,disphy=serialize(disph_data[:,0],disph_data[:,1])
chdfx,chdfy=serialize(chdf_data[:,0],chdf_data[:,1])

#ch_func=ip.interp1d(np.log10(chdfx),np.log10(chdfy),kind='slinear')
#disph_func=ip.interp1d(np.log10(disqkx),np.log10(disqky),kind='slinear')
#disqk_func=ip.interp1d(np.log10(disphx),np.log10(disphy),kind='slinear')
ch_func=ip.interp1d(np.log10(chdfx),np.log10(chdfy),kind='slinear')
disph_func=ip.interp1d(np.log10(disphx),np.log10(disphy),kind='slinear')
disqk_func=ip.interp1d(np.log10(disqkx),np.log10(disqky),kind='slinear')
en_arr=np.logspace(2,8,1000)

disxs=(10**(disph_func(4.))+10**(disqk_func(4.)))*1e4*1e36

ph_xs=10**(disph_func(4.))*1e4*1e36
qk_xs=10**(disqk_func(4.))*1e4*1e36
print ("DIS cross-section for 10TeV neutrino is (pb): ", disxs)

print ("Photon: ", ph_xs)
print ("Quark: ", qk_xs)
sumxs=10**(ch_func(np.log10(en_arr)))+10**(disph_func(np.log10(en_arr)))+10**(disqk_func(np.log10(en_arr)))

plt.figure()
#plt.plot(totalx,totaly,lw=1,label='Total XS')
#plt.scatter(chdfx,chdfy,s=5,label='CHDF data')
#plt.scatter(disqkx,disqky,s=5,label='DIS data')
#plt.scatter(disphx,disphy,s=5,label='DIS data')

plt.plot(en_arr,10**(ch_func(np.log10(en_arr))),lw=1,label='Coherent + Difffractive')
plt.plot(en_arr,10**disph_func(np.log10(en_arr)),lw=1,label='Photon Initiated DIS')
plt.plot(en_arr,10**disqk_func(np.log10(en_arr)),lw=1,label='Quark Initiated DIS')

plt.plot(en_arr,sumxs,lw=1,label='Total Cross-section')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\sigma/E_{\nu}$ $[cm^2GeV^{-1}]$')
plt.xlabel('Nuetrino Energy (GeV)')
plt.title(r'Neutrino trident cross-section for $\nu_{\mu}\rightarrow \nu_{\mu}+\mu^-+\mu^+$ on $H_2O$')
plt.legend(fontsize=5)
plt.grid(which='both',alpha=0.3,ls='--')
plt.savefig(pltloc+'xs_check.pdf')
