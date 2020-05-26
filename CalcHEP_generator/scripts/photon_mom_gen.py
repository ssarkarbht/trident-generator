#!/bin/python

import numpy as np
import argparse

parser = argparse.ArgumentParser("Neutrino energy")
parser.add_argument('--energy', '-e', type=float)
args=parser.parse_args()

nu_en = args.energy
#print (args.energy)

#following limit is for tridents
m2jk = 0.044944
qmin = m2jk/4.0/nu_en
qmax = 1.
qlogmin = np.log10(qmin)
#qlogmax = np.log10(1.)
qlogmax = np.log10(qmax)
x = np.logspace(qlogmin,qlogmax,30)


#W-Boson resonance energy range:
cmres=6.5e3
cmresmin=3e3
cmresmax=9e3
qresmin = cmresmin/4.0/nu_en
qresmax = cmresmax/4.0/nu_en
lqresmin=np.log10(qresmin)
lqresmax=np.log10(qresmax)
y=np.logspace(lqresmin,lqresmax,20)

xp=y[np.where(y<qmax)]

x=sorted(np.append(x,xp))

for i in x:
        print ("%.4e" % i)
