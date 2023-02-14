#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:33:59 2023

Creates a Probability density function of a ligand to each residue.
Useful for estimating contacts between two molecules or domains.

Checked code with VMD (Note VMD uses center of geometry, not center of mass)
and results match.

@author: bab226
Version: 1
"""

import MDAnalysis as mda
from MDAnalysis.analysis import contacts

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

top = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/cluster/cluster_helix/final_edit.gro"
#xtc = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/cluster/cluster_helix/clusters.pdb001.pdb" # Test data
xtc = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/final_100-200ns_edit.xtc"
u = mda.Universe(top, xtc)

contacts_per_res = []
residues = np.array(list(range(1,38)))
    
for resnumber in range(1, 38):
    res_contacts = 0
    
    for ts in u.trajectory[:]:
        
        sel_res = "(resnum %s)" %(resnumber)
        sel_adm = "(index 537-645)"
        
        res = u.select_atoms(sel_res)
        adm = u.select_atoms(sel_adm)
        
        dist = res.center_of_mass() - adm.center_of_mass()
        dist = np.linalg.norm(dist)   #Normalize vector
    
        if dist < 10:
            res_contacts = res_contacts + 1
            
        else: 
            res_contacts = res_contacts
        
    contacts_per_res.append(res_contacts)

contact_density = contacts_per_res/np.sum(contacts_per_res)
#%%
np.savez_compressed('/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/contacts/contact-density.npz', a=residues, b=contact_density)
#%%

plt.plot(residues, contact_density)
plt.ylabel('# Contacts')
plt.xlabel("Resid #")
