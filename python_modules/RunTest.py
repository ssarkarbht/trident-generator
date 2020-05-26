#!/bin/python3

import numpy as np


print ("Choose options for testing modules: [1:Muon Range, 2:Read LHE files]")
op=input("Input test number: ")
if op=='1':
	print ("Testing muon range module...")
	import muon_range as mr

#print ("testing energy: 2e2 GeV")
	var=mr.range_func(2e2)
	if np.isnan(var):
		print ("Muon range module FAILED")
	else:
		print ("Muon range module PASSED")

elif op=='2':
	print ("Testing Read_LHE module...")
	import Read_LHE as rl
	from glob import glob
	files=glob("*.lhe")
	if len(files)==0:
		print ("Testing FAILED: No LHE file present in the current directory!")
	else:
		evt=rl.LHEEVE(files[0])
		mum=evt.mm_mom()
		if len(mum):
			print ("Read LHE file module PASSED!")
		else:
			print ("Read LHE file module FAILED!")

