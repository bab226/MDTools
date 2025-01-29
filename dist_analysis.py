#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:33:59 2023

Determines the closest residue to a bound ligand.
Useful for estimating determining ligand affinity to a target.

Checked code with VMD (Note VMD uses center of geometry, not center of mass)
and results match.

@author: bab226
Version: 1
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
import numpy
import os
import sys

path = "./analysis"
isExist = os.path.exists(path)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(path)
    print("Analysis directory is created!")
else:
    print("Analysis Directory Already Exists!")
    
target = sys.argv[1]
inhibitor = sys.argv[2]
replicate = sys.argv[3]
# target = "hiapp_cterm_beta"
# inhibitor = "xnfxvh"
# replicate = "r1"

top = "./complexes/%s/%s/%s/" %(target, inhibitor, replicate) + "nvt_noj.gro"
#top = "./testME/nvt_noj.gro"
xtc = "./complexes/%s/%s/%s/" %(target, inhibitor, replicate) + "nvt_noj_res_whole.xtc"
#xtc = "./testME/nvt_noj_res_whole.xtc"
u = mda.Universe(top, xtc)

frames = len(u.trajectory)

contacts_per_res = []
#residues = np.array(list(range(1,20)))

avg_dist_arr = []
df_full_dist = pd.DataFrame()

#sd_dist_arr = []

if "cterm" in target:  
    res_start = 18   # C-term starts at residue 19
    c = 0   # C-term starts at ACE-0
else:
    res_start = 1
    c = 1   # N-term starts at LYS-1 (therefore add 1 to resnumber)

if "adm" in inhibitor:
    residue_num = len(u.residues)-1     # Inhibitor is last residue.
else:
    residue_num = len(u.residues)-6     # Inhibitors are last six residues.

#print(residue_num)
resid = list(range(res_start, res_start + residue_num))
  
resnumber = 10   # Picking one residue in target arbitrarily

#res_contacts = 0

residue = str(u.residues[resnumber])[9:12]

### Select residues from GRO file based on residue and confirmed by index number.
### Errors might occur here so check closely with GRO file!

if "hiapp_cterm" in target:  
    sel_res = "(resnum %s and index 0-266)" %(resnumber+c) #0-266 is the indices for protein in Human C-term.

elif "hiapp_nterm" in target:  
    sel_res = "(resnum %s and index 0-331)" %(resnumber+c) #0-332 is the indices for protein in Human N-term.

elif "riapp_nterm" in target:  
    sel_res = "(resnum %s and index 0-338)" %(resnumber+c) #0-339 is the indices for protein in Rat N-term.
    
elif "riapp_cterm" in target:  
    sel_res = "(resnum %s and index 0-272)" %(resnumber+c) #0-273 is the indices for protein in Rat C-term.

print(sel_res)

# Get initial coordinates of inhibitor 

for ts in u.trajectory[0]:
    if "adm" in inhibitor:
        if "cterm" in target:  
            inh = u.residues.indices[20]  # ADM in C-term system is 20
        else:
            inh = u.residues.indices[23]  # ADM in N-term system is resid 24
            
        start = np.min(inh)
        end = np.max(inh)
        sel_inh = "index %s-%s" %(start, end)
    else:
        if "cterm" in target:  
            inh = u.residues.indices[20:26]  # Hexapeptide inhibitor in N-term is resid 20:26 (both C-term GRO and python starts at 0)
        else:
            inh = u.residues.indices[23:29]  # Hexapeptide inhibitor in N-term system is resid 24:30 (N-term GRO file starts at 1, but python array starts at 0)
            
        start = np.min(inh[0])
        end = np.max(inh[5])
        sel_inh = "index %s-%s" %(start, end)
    
    #print(sel_inh)
    res = u.select_atoms(sel_res)
    
    inh = u.select_atoms(sel_inh)
    
    origin = [50, 50, 50]

#Itarate through frames from simulation and compare distance from COM of origin to inhibitor COM over all frames.

# Define start and end for binning data
frame_arr = []
dist_sq_arr = []
residue_arr = []

# Go through all frames iteratively
for ts in u.trajectory[0:frames]:
    
    if "adm" in inhibitor:
        if "cterm" in target:  
            inh = u.residues.indices[20]  # ADM in C-term system is 20
        else:
            inh = u.residues.indices[23]  # ADM in N-term system is resid 24
            
        start = np.min(inh)
        end = np.max(inh)
        sel_inh = "index %s-%s" %(start, end)
    else:
        if "cterm" in target:  
            inh = u.residues.indices[20:26]  # Hexapeptide inhibitor in N-term is resid 20:26 (both C-term GRO and python starts at 0)
        else:
            inh = u.residues.indices[23:29]  # Hexapeptide inhibitor in N-term system is resid 24:30 (N-term GRO file starts at 1, but python array starts at 0)
            
        start = np.min(inh[0])
        end = np.max(inh[5])
        sel_inh = "index %s-%s" %(start, end)
    
    #print(sel_inh)
    res = u.select_atoms(sel_res)
    inh = u.select_atoms(sel_inh)
    
    dist = origin - inh.center_of_mass()   # Displacement
    dist = np.linalg.norm(dist)
    dist_sq = dist**2   #Displacement-squared
    
    dist_sq_arr.append(dist_sq)
    frame_arr.append(ts.frame)
    
print(res[-1:])
print(inh[-1:])

time_arr = np.array(frame_arr) * ts.dt #Multiply frames by timestep
residue_arr.append(residue)
    # print(residue)
    # print(avg_dist_arr[resnumber])
    
cols = {'time': time_arr,
        'dist': dist_sq_arr,
        #'resid' : resid,
        }

df_dist = pd.DataFrame(cols)

path = "./analysis/distance"
isExist = os.path.exists(path)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(path)
    print("Distance directory is created!")
else:
    print("Distance Directory Already Exists!")
       
   
df_dist.to_csv("./analysis/distance/%s-%s-%s-dist-time.csv" %(target, inhibitor, replicate))

         


