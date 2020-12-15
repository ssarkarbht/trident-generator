#!/bin/python

'''
Author: Sourav Sarkar
Date: October 8, 2020
Objective: This scritps injects muon neutrino around the detector
	with ranged mode and creates CC muons in the interaction
'''

from I3Tray import *
from icecube import icetray, dataclasses, dataio
from icecube import phys_services
from icecube import LeptonInjector
from icecube.icetray import I3Units
from icecube.hdfwriter import I3HDFWriter

seed=29

#creating icetray
tray=I3Tray()

# Add random service, earth model and source of blank DAQ frames
randomService=phys_services.I3GSLRandomService(seed=seed)
tray.context["I3RandomService"] = randomService
tray.AddService("I3EarthModelServiceFactory", "Earth")
tray.AddModule("I3InfiniteSource", "TheSource", Stream=icetray.I3Frame.DAQ)


#Path to differential cross-section data: Using standard cross-section files from NuGen
xs_folder = "/cvmfs/icecube.opensciencegrid.org/data/neutrino-generator/cross_section_data/csms_differential_v1.0/"

tray.AddModule("RangedLeptonInjector",
	DoublyDifferentialCrossSectionFile=xs_folder+'/dsdxdy_nu_CC_iso.fits',
	TotalCrossSectionFile             =xs_folder+"/sigma_nu_CC_iso.fits",
	FinalType1                        =dataclasses.I3Particle.ParticleType.MuMinus,
	FinalType2                        =dataclasses.I3Particle.ParticleType.Hadrons,
	EarthModel                        ="Earth",
	EndcapLength                      =1200*I3Units.meter,
	InjectionRadius                   =1200*I3Units.meter,
#	CylinderRadius                    =800.*I3Units.meter,
#	CylinderHeight                    =800.*I3Units.meter,
	MinimumZenith                     =90 * I3Units.deg,
	MaximumZenith                     =180 * I3Units.deg,
	PowerLawIndex                     = 1.,
	MaximumEnergy                     =(1e8)*I3Units.GeV,
	MinimumEnergy                     =100.*I3Units.GeV,

	NEvents                           =20000,
	RandomService                     ="I3RandomService")

# a little function to add the event number and seed to each frame 
event_id = 1
def get_header(frame):
    global event_id
    header          = dataclasses.I3EventHeader()
    header.event_id = event_id
    header.run_id   = seed
    frame["I3EventHeader"] = header

    event_id += 1

tray.AddModule(get_header, streams = [icetray.I3Frame.DAQ])
# write two files: one the configuration file for use with LeptonWeighter
tray.AddModule("InjectionConfigSerializer", OutputPath ="dataset_03/numu_batch5_config.lic")

# Write the datafile 
tray.AddModule("I3Writer", Filename="dataset_03/numu_batch5.i3",
        Streams=[icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation,
                 icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

tray.Execute()
tray.Finish()

