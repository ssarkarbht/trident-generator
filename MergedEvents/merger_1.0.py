#!/bin/python
from icecube import dataio,dataclasses
import numpy as np
from glob import glob
import Read_LHE as rl

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-b", "--BatchNumber", dest="BATCH", type=int)

(options, args) = parser.parse_args()
batchnum=options.BATCH


input_dict={1:'batch_1/',2:'batch_2/',3:'batch_3/',4:'batch_4/',5:'batch_5/'}
output_dict={1:'i3input_D01_B1.txt',2:'i3input_D01_B2.txt',3:'i3input_D01_B3.txt',4:'i3input_D01_B4.txt',5:'i3input_D01_B5.txt'}

i3fileloc='/home/ssarkar/NTPGenerator/LeptonInjector/'

i3file=dataio.I3File(i3fileloc+'numu_batch2.i3')

eventfileloc='/data/icecube/ssarkar/dataio/numu_dataset_02/'+input_dict[batchnum]

datafiles=sorted(glob(eventfileloc+'event_id*'))

def fill_event(mctree, eventfile, etype):
	if etype=='lhe':
		geo=np.array([-1.])
		nu_i=np.array([1.,14.])
		nu_f=np.array([2.,14.])
		mm_f=np.array([2.,13.])
		mp_f=np.array([2.,-13.])
		pp_f=np.array([2.,2212.])

		evt_pos=np.array([mctree[0].pos.x,mctree[0].pos.y,mctree[0].pos.z])
		evt_dir=np.array([mctree[0].dir.theta,mctree[0].dir.phi])
		geo=np.append((geo),evt_pos)
		geo=np.append(geo,evt_dir)

		nu_mom=eventfile.primary_nu()
		nu_i=np.append(nu_i,nu_mom)

		nu_mom1=eventfile.nm_mom()
		nu_f=np.append(nu_f,nu_mom1)

		mm=eventfile.mm_mom()
		mm_f=np.append(mm_f,mm)

		mp=eventfile.mp_mom()
		mp_f=np.append(mp_f,mp)

		pp=eventfile.add_recoil()
		pp_f=np.append(pp_f,pp)

		event=np.vstack((geo,nu_i,nu_f,mm_f,mp_f,pp_f,np.zeros(6)))

		return event

	if etype=='hepmc':
		geo=np.array([-1.])
		evt_pos=np.array([mctree[0].pos.x,mctree[0].pos.y,mctree[0].pos.z])
		evt_dir=np.array([mctree[0].dir.theta,mctree[0].dir.phi])
		geo=np.append((geo),evt_pos)
		geo=np.append(geo,evt_dir)

		nu_i=eventfile.inu_mom()
		nu_f=eventfile.onu_mom()
		mm_f=eventfile.mm_mom()
		mp_f=eventfile.mp_mom()
		shower=eventfile.shower_mom()

		event=np.vstack((geo,nu_i,nu_f,mm_f,mp_f,shower,np.zeros(6)))

	return event

event_array=np.zeros(6)
print ("Starting event processing...")
fileindex=0
while i3file.more():
	frame=i3file.pop_daq()
	evtid=frame['I3EventHeader'].event_id
	mctree=frame['I3MCTree']
	df=datafiles[fileindex]
	if df[-5:]=='hepmc':
		etype='hepmc'
		eid=df[-11:-6]
		hep=rl.HEPMC(df)
		if hep.evtnum==0:
			print ("Skipping empty HepMC file: ", df)
			fileindex+=1
			continue
		event=fill_event(mctree,hep,etype)
		event_array=np.vstack((event_array,event))
		fileindex+=1
	else:
		etype='lhe'
		eid=df[-12:-7]
		index=np.random.randint(1000)
		lhe=rl.LHEEVE(df,index)
		event=fill_event(mctree,lhe,etype)
		event_array=np.vstack((event_array,event))
		fileindex+=1
	if evtid%1000==0:
		print ("No. of events processed and matched: ", evtid, eid)

np.savetxt('/home/ssarkar/NTPGenerator/Event_Merger/dataset_02/'+output_dict[batchnum],event_array)
