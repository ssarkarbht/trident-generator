#!/bin/python3

'''
Author: Sourav Sarkar
Date: 22 April, 2020
Description: This script reads the lhe event file generated in CalcHEP
		and converts the event properties in a hdf5 file
'''

import numpy as np
import h5py
from glob import glob
import Read_LHE as rl
import argparse

parser = argparse.ArgumentParser("Neutrino energy")
parser.add_argument('--energy', '-e', type=str)
args=parser.parse_args()

nu_en = args.energy

loc='/data/user/ssarkar/CalcHEP/event_database_chdf'

#For testing, only calling one neutrino energy directory

enudir="En_"+nu_en
datadir='eventfiles'
if np.float_(nu_en)>=1e5:
	datadir='eventfiles_add'

infiles=glob(loc+'/'+enudir+'/'+datadir+'/*.lhe')

h5f=h5py.File(enudir,"w")

for i in infiles:
	q=i[-14:-4]
	evt=rl.LHEEVE(i)
	mum=evt.mm_mom()
	mup=evt.mp_mom()
	opang=evt.op_ang(mum,mup)[:,np.newaxis]
	p1=mum[:,3][:,np.newaxis]
	p2=mup[:,3][:,np.newaxis]
	dataset=np.hstack((p1,p2,opang))
	h5f.create_dataset(q,data=dataset)
#print (h5f.keys())

