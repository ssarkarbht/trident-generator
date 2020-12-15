from I3Tray import *
from os.path import expandvars
import os, sys
from optparse import OptionParser

from icecube import icetray, dataclasses, dataio, phys_services
from icecube import clsim

parser=OptionParser()
parser.add_option("-g", "--InputGCD", dest="GCD", type=str)
parser.add_option("-i", "--InputFile", dest="INPUT", type=str)
parser.add_option("-o", "--OutputFile", dest="OUTPUT", type=str)
parser.add_option("-s", "--RandomSeed", dest="SEED", type=int)
(options, args) = parser.parse_args()

I3Seed = 12345+options.SEED

randomService = phys_services.I3SPRNGRandomService(
        seed = I3Seed,
        nstreams = 1000,
        streamnum = 1)


#filter module (Point of Closest Approach cut)
def run_clsim(frame):
	return frame['PCA_cut']==icetray.I3Bool(True)


tray = I3Tray()

tray.AddModule("I3Reader", "reader", FilenameList = [options.INPUT])
#tray.AddModule(run_clsim,"filter",streams=[icetray.I3Frame.DAQ])

tray.AddSegment(clsim.I3CLSimMakeHits, "makeCLSimHits",
		DOMOversizeFactor=1,
		MCTreeName="I3MCTree_postMuonProp",
#		MMCTrackListName="MMCTrackList",
		If=run_clsim,
		GCDFile=options.GCD,
		RandomService=randomService,
		IgnoreSubdetectors=["IceTop"],
		UseCPUs=True)#,
#		UseGPUs=True)

tray.AddModule("I3Writer", "writer", Filename = options.OUTPUT)

tray.AddModule("TrashCan", "trash")

tray.Execute(1)
tray.Finish()

