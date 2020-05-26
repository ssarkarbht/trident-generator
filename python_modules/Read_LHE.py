#!/bin/python
'''
Author: Sourav Sarkar
Date: 21-01-2020
Objective: Module to call event property functions that read the LHE files
		and returns the values of each event

WARNING!! Only SPACE indentation is used in this script, DO NOT use TAB indentation
Mixing of TAB and SPACE indentation does not go well in python scripts, apologize for the inconvenience
'''
import numpy as np
from scipy.interpolate import interp1d
import muon_range as mr

class LHEEVE:
        #Initializing the class with loading the datafile
        def __init__(self,fname):
                block = '<event>'
                f_in = open(fname,'r')
                lines = np.asarray(f_in.readlines())
                self.mlines = np.char.strip(lines)
                self.event_index = np.where(self.mlines==block)[0]
        #Calculatating Q^2 from incoming and outgoing nucleus (quarks)
        def Q2(self):
                i_pos = self.event_index+3 #incoming nucleus
                f_pos = self.event_index+7 #outgoing nucleus

                #data saved as array of lists
                p_i = np.char.split(self.mlines[i_pos])
                p_f = np.char.split(self.mlines[f_pos])
                #convert it into list of lists
                p_ii = p_i.tolist()
                p_ff = p_f.tolist()
                #convert it back to array of arrays (and keeping only 4-momentum)
                p_i = np.asarray(p_ii).astype(np.float)[:,6:10]
                p_f = np.asarray(p_ff).astype(np.float)[:,6:10]

                #Compute the Q2 values on the sorted arrays
                qx = p_i[:,0]-p_f[:,0]
                qy = p_i[:,1]-p_f[:,1]
                qz = p_i[:,2]-p_f[:,2]
                qt = p_i[:,3]-p_f[:,3]

                q2arr = -1.*(qt*qt-qx*qx-qy*qy-qz*qz)
                #now return the 1-dim Q^2=-q^2 array
                return q2arr

        #Calculate rest mass of any particle
        def mass(self, mom4):
                dot=mom4[:,3]*mom4[:,3]-mom4[:,0]*mom4[:,0]-mom4[:,1]*mom4[:,1]-mom4[:,2]*mom4[:,2]
                return (np.sqrt(dot))
        #extract outgoing neutrino 4-momentum
        def nm_mom(self):
                nm_pos = self.event_index+4 #outgoing nu_mu
                nm_p = np.char.split(self.mlines[nm_pos])
                nm_pp = nm_p.tolist()
                nm_p = np.asarray(nm_pp).astype(np.float)[:,6:10]
                return nm_p
        #extract outgoing mu- 4-momentum
        def mm_mom(self):
                mm_pos = self.event_index+5 #outgoing mu-
                mm_p = np.char.split(self.mlines[mm_pos])
                mm_pp = mm_p.tolist()
                mm_p = np.asarray(mm_pp).astype(np.float)[:,6:10]
                return mm_p
        #extract outgoing mu+ 4-momentum
        def mp_mom(self):
                mp_pos = self.event_index+6 #outgoing mu+
                mp_p = np.char.split(self.mlines[mp_pos])
                mp_pp = mp_p.tolist()
                mp_p = np.asarray(mp_pp).astype(np.float)[:,6:10]
                return mp_p

        #Calculate the opening angle between two outgoing muon tracks
        def op_ang(self, mm, mp):
                #dot product between the two vectors
                dot = mm[:,0]*mp[:,0]+mm[:,1]*mp[:,1]+mm[:,2]*mp[:,2]
                norm_m = np.linalg.norm(mm[:,0:3],axis=1)
                norm_p = np.linalg.norm(mp[:,0:3],axis=1)
                ct = dot/norm_m/norm_p
# Some cos(theta) values go out of bound due to rounding error
                outidx=np.where(abs(ct)>1.)
                ct[outidx]=1.
                theta = np.arccos(ct)
                return theta #in radian
        #Calculate the invariant mass of two outgoing particles 4-momentum
        #in case of W-boson invariant mass: outgoing numu & mu+
        def inv_mass(self, nm, mp):
                mu_mass = 0.106 #muon mass in GeV
                #dot product of 3-momentum
                dot = mp[:,0]*nm[:,0]+mp[:,1]*nm[:,1]+mp[:,2]*nm[:,2]
                invm = mu_mass**2+2.*(mp[:,3]*nm[:,3]-dot)
                return invm
        #Calculate the minimum track length of the two muon track
        def min_track(self, mm, mp):
                min_en = np.minimum(mm[:,3],mp[:,3])
                #call the muon stopping power function
                return mr.range_func(min_en)
