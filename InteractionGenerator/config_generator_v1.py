#!usr/bin/env python

'''
Author: Sourav Sarkar
Date: Nov 5, 2020
Description: This script takes the LI I3 file and generates the event generation
		configuration file for running CalcHEP and MadGraph
	Configuration generator samples necessary interaction type, nucleus type, photon energy
	and stores into output HDF5 file. For backtracking purpose, event id from I3 file is 
	also copied to HDF5 file for later use of pulling additional information from I3files
	and creating final Trident events with full event info (complete I3MCTree and 
	I3MCWeightDictionary)

'''
import numpy as np
from icecube import dataclasses,dataio,LeptonInjector,phys_services
import generator as gn
from optparse import OptionParser
import h5py as h5
from numpy.random import default_rng


#File location for the Lepton Injector i3 files
LI_filepath='/data/user/ssarkar/TridentProduction/simulation/LeptonInjector/dataset_03/'
out_filepath='dataset_03/'
#Add command-line argument options
parser=OptionParser()
parser.add_option("-s", "--RandomSeed", dest="SEED", type=int, default=12345)
parser.add_option("-n", "--EventNumber", dest="NUM", type=int, default=10000)
parser.add_option("-i", "--Inputfilename", dest="INFILE", type=str, default='sample_out.i3')
parser.add_option("-o", "--OutputfileName", dest="OUTFILE",type=str, default='sample_output.h5')


#parse the parameters
(options, args)=parser.parse_args()

#load the Trident cross-section files
#file position: 0: EPA xs, 1: coherent+diffractive, 2: DIS_quark, 3: DIS_photon
xsfiles=['cross-section_data/EPA_CM_SM_XS_DATA.txt', 'cross-section_data/new_beacom_chdf.csv',
        'cross-section_data/new_beacom_DIS_quark.csv', 'cross-section_data/new_beacom_DIS_photon.csv']

#prepare the Oxygen and Hydrogen (coherent+diffractive) weight factors
log_wf=gn.weight_function(['OX_weight_factor_v1.txt','H_weight_factor_v1.txt'])

#call the numpy random generator
seed=options.SEED+12345
#numpy random generator
np_rng=default_rng(seed)
#icecube random generator
ic_rng=phys_services.I3GSLRandomService(seed=seed)

#get the other parameters
evtnum=options.NUM
outfile=out_filepath+options.OUTFILE
infile=LI_filepath+options.INFILE


#Create input and output file object
i3file=dataio.I3File(infile,'r')
h5file=h5.File(outfile,'w')

#Two branches for two different interaction regimes 
chdf_branch=h5file.create_group("CH+DF_interactions")
dis_branch=h5file.create_group("DIS_interactions")

#structured array to fill the event config parameters with (for two different regimes)
chdf_arr=np.zeros(evtnum,dtype=[('event_id',int),('nu_pdgid',int),('nu_energy',float),('photon_energy',float)])
dis_arr=np.zeros(evtnum,dtype=[('event_id',int),('nu_pdgid',int),('nu_energy',float),('nucleus_type', int),('DIS_type',int)])

#call the generator class object
config=gn.generator(np_rng,xsfiles)

#array iterators for both containers
chdf_iter=0
dis_iter=0

#counter for logging purpose
#counter=1

#Loop throught the events in I3Files
while i3file.more():
#	if counter%1000==0:
#		print ("Number of events processed: ", counter)
#	counter+=1

#extrtact information from i3 file
	frame=i3file.pop_daq()
	evtid=frame["I3EventHeader"].event_id
	if evtid%1000==0:
		print ("Number of events processed: ", evtid)
	mctree=frame['I3MCTree']
	prim=mctree[0]
	pdg=prim.pdg_encoding
	en=prim.energy
#sample interaction type
	int_type=config.sample_interaction(en)
#Fill configuration parameters based on the interaction type
	if int_type==1:
		photon=gn.EPAPhoton(en,ic_rng,log_wf)
		phen=photon.sampling(1)
		param=(evtid,pdg,en,phen)
		chdf_arr[chdf_iter]=param
		chdf_iter+=1
	elif int_type==2:
		dis_type=config.sample_DIS(en)
		n_type=config.sample_nucleus()
		param=(evtid,pdg,en,n_type,dis_type)
		dis_arr[dis_iter]=param
		dis_iter+=1
# If interaction type is different from the above two, something is wrong!
	else:
		print ("Invalid Interaction Type!! Check for bugs! ", int_type)

		

#throw away the empty rows from each branch containers
chdf_cut=np.delete(chdf_arr,np.s_[chdf_iter:])
dis_cut=np.delete(dis_arr,np.s_[dis_iter:])

#finally store the arrays into hdf5 file
chdf_branch.create_dataset("interaction_param",data=chdf_cut)
dis_branch.create_dataset("interaction_param", data=dis_cut)

