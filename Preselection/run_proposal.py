#!/usr/bin/env python

from icecube import dataclasses, dataio, phys_services, icetray, simclasses,sim_services

from I3Tray import *

from icecube.simprod.segments import PropagateMuons
from optparse import OptionParser
import numpy as np

parser=OptionParser()
parser.add_option("-i", "--InputFile", dest="INPUT", type=str)
parser.add_option("-s", "--RandomSeed", dest="SEED", type=int)
parser.add_option("-o", "--OutputFile", dest="OUTPUT", type=str)

(options, args) = parser.parse_args()


#prepare the poca cut function
center=dataclasses.I3Position(0.,0.,0.)
calc=phys_services.I3Calculator

fileloc='/data/user/ssarkar/TridentProduction/simulation/MergedEvents/dataset01/'
infile=fileloc+options.INPUT
seed=options.SEED
outfile=options.OUTPUT


def proj_poca_distance(pos,cutvalue):
	xy=np.sqrt(pos.x**2+pos.y**2)
	z=abs(pos.z)
	if xy<=cutvalue and z<=cutvalue:
		return True
	return False

def put_pocacut(frame):
	global center
	cutvalue=500.

	mctree=frame['I3MCTree_postMuonProp']
	prim=mctree.primaries[0]
	seco=mctree.get_daughters(prim)

	mm=seco[1]
	mp=seco[2]

	mmpocapos=calc.closest_approach_position(mm,center)
	mppocapos=calc.closest_approach_position(mp,center)

	if proj_poca_distance(mmpocapos,cutvalue) and proj_poca_distance(mppocapos,cutvalue):
		frame['PCA_cut']=icetray.I3Bool(True)
	else:
		frame['PCA_cut']=icetray.I3Bool(False)

tray=I3Tray()

randomService = phys_services.I3SPRNGRandomService(
		seed=seed,
		nstreams=20000,
		streamnum=1)

tray.AddModule("I3Reader", 'reader', FileNameList= ['GeoCalibDetectorStatus_2016.57531_V0.i3.gz', infile])

tray.AddSegment(PropagateMuons, 'MuonPropagation',
		RandomService=randomService,
		InputMCTreeName='I3MCTree',
		OutputMCTreeName='I3MCTree_postMuonProp')
tray.AddModule(put_pocacut,'pocacut',Streams=[icetray.I3Frame.DAQ])
tray.AddModule("I3Writer", 'writer', FileName=outfile)
tray.AddModule("TrashCan", "simplytrash")

tray.Execute()
tray.Finish()


