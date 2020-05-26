This directory contains generic python modules that are used across different
tools and purposes:

----------------------------------------------
Module: Muon track range in Ice

PythonFile:muon_range.py

Datafile: mue_water_ice.txt

DataSource: PDG particle data

HighEnergy Extrapolation: PDG data has the highest Muon range of 1e6 GeV
However, our muon track energy can reach as high as 1e8 GeV. In order to make a rough estimation of the muon track length in ice (rough because at high energies muon mainly loose energy via stochastic loss which has heavy statistical fluctuation)
----------------------------------------------
----------------------------------------------
Module: Read LHE eventfiles

PythonFile:Read_LHE.py

Input: LHE eventfile, in order to test this module, a lhe file has top be present in the same directory
----------------------------------------------

