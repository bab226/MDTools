#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:35:27 2022

Helper functions for MDTraj Analysis

quant_dssp and bootstrapping in dssp_w_error done by 
Dr. Gregory 'Beggory-Sherlock' Gomes

Adapted by
@author: bab226
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import math
import pandas as pd
import mdtraj as md

from random import seed
from random import randint
from random import random
from random import uniform

SMALLEST_SIZE = 2.0 #was 1.5
SMALL_SIZE = 12/2.25
MEDIUM_SIZE = 18/2.25
BIGGER_SIZE = 20/2.25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

path = '/Users/bab226/Documents/yale_research/iapp/md_sims/iapp_adm116/new_files/analysis/'
aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def quant_dssp(traj,w=np.array([])):
    '''
    From Dr. Gregory 'Beggory-Sherlock' Gomes
    
    Parameters
    ----------
    traj : mdtraj trajectory
    w : weights, optional
        DESCRIPTION. The default is np.array([]).

    Returns
    -------
    coil : count of coil
    helix : count of helix
    strand : count of strand

    '''
    coil=np.zeros(traj.topology.n_residues)
    helix=np.zeros(traj.topology.n_residues)
    strand=np.zeros(traj.topology.n_residues)
    
    #Normalization to 1
    if not w.any():
        w=np.ones(len(traj))
        w=w/len(w)
    
    for i in range(len(traj)):
        dssp=md.compute_dssp(traj[i])[0]
        for j in range(len(dssp)):
            if 'C' in dssp[j]:
                coil[j]+=w[i]
            elif 'H' in dssp[j]:
                helix[j]+=w[i]
            else:
                strand[j]+=w[i]
    return coil, helix, strand


def dssp_w_error(traj,w=np.array([])):
    '''
    From Dr. Gregory 'Beggory-Sherlock' Gomes
    
    Parameters
    ----------
    traj : mdtraj trajectory
    w : weights, optional
       DESCRIPTION. The default is np.array([]).
       
    Returns
    -------
    C : coil
    C_std : coil sd
    S : strand
    S_std : strand sd
    H : helix
    H_std : helix sd
 
    '''
	#Use bootstrapping to estimate the statistical error in fraction coil, helix, strand

	#If you don't include a vector of weights, it will give each conformation an equal weight (e.g., 1/Nconformation)
    #This is what you want for standard MD
    
    if not w.any():
        w=np.ones(len(traj))
        w=w/len(w)    

    num_straps = 1 #NUMBER OF BOOTSTRAPS TO DO. CHOOSE LOW NUMBER FOR TESTING (e.g., 5) CHOOSE HIGH NUMBER FOR PROPER ESTIMATE (e.g., 1000)
    C_array = []
    H_array = []
    S_array = []
    for i in range(num_straps):
        print('Bootstrap %d of %d' %(i,num_straps) )
        wsi = np.random.choice(len(w), size=len(w), p=w/np.sum(w), replace=True) #Doing the bootstrapping! Choose N samples with replacement with probability w (does weighting and bootstrapping)
        C_0, H_0, S_0 = quant_dssp(traj[wsi])
        C_array.append(C_0)
        H_array.append(H_0) 
        S_array.append(S_0)        
        
    C_array = np.array(C_array) 
    H_array = np.array(H_array)
    S_array = np.array(S_array)    


    C = np.mean(C_array,axis = 0)
    C_std = np.std(C_array, axis = 0, ddof = 1)
    S = np.mean(S_array,axis = 0)
    S_std = np.std(S_array, axis = 0, ddof = 1)
    H = np.mean(H_array,axis = 0)
    H_std = np.std(H_array, axis = 0, ddof = 1)

    return C, C_std, S, S_std, H, H_std



def dssp_plot(traj, H, H_std, S, S_std, H2, H2_std, S2, S2_std, path, title, labels, ss):
    '''
    
    Parameters
    ----------
    aa : amino acids
    C : coil
    C_std : coil sd
    S : strand
    S_std : strand sd
    H : helix     
    H_std : helix sd
    title : string
    labels : array of strings
    ss : string  ## 'H', 'S', or 'all'

    Returns
    -------
    None.

    '''
    aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
   
    aa=[x.resSeq for x in traj.topology.residues]
    # print(aa)
    aa_name = [aa_dict[x.name] for x in traj.topology.residues]
    print(aa_name)
    n = ""
    n = [n.join(aa_dict[x.name]) for x in traj.topology.residues]
    print("".join(str(x) for x in n))

    psites = np.array([])  #Will put vertical lines at these sites to emphasize the residues.
    psites = psites.astype(int)    
    
    wth = 2.3 *1.2 #originall *1
    hgt = wth/1.25    

    lw = 1
    ms = 1 

    psites = np.array([]) - 1 #Will put vertical lines at these sites. In IAPP these are cysteine disulfide bond sites.
    psites = psites.astype(int)    

    # plt.figure()
    # fig, ax = plt.subplots(figsize=(10, 8))
    # plt.subplots_adjust(left=0.15, right=0.95, top=0.83, bottom=0.12)
    
    fig, ax = plt.subplots(figsize=(wth, hgt))
    plt.subplots_adjust(left=0.2, right=0.85, top=0.83, bottom=0.18)
    # plt.subplots_adjust(left=0.2, right=0.9, top=0.83, bottom=0.18)

    for p in psites:
		#put vertical lines at defined sites
        ax.axvline(np.array(aa)[p], color = 'k', alpha = 0.2, lw = lw)
        
    label1 = labels[0]
    label2 = labels[1]
    
	#The actual error bar plots
    if ss == "H":
        label = "% Helix"
        ax.errorbar(aa, H, yerr = H_std, color = 'orange', linestyle='--', label = label1, lw = lw, ms = ms)
        ax.errorbar(aa, H2, yerr = H2_std, color = 'orange', linestyle='-', label = label2, lw = lw, ms = ms)
    elif ss == "S":
        label = "% Sheet"
        ax.errorbar(aa, S, yerr = S_std, color = 'green', linestyle='--', label = label1, lw = lw, ms = ms)
        ax.errorbar(aa, S2, yerr = S2_std, color = 'green', linestyle='-', label = label2, lw = lw, ms = ms)
    else:
        label = "SSP"
        ax.errorbar(aa, H, yerr = H_std, color = 'orange', linestyle='--', label = label1, lw = lw, ms = ms)
        ax.errorbar(aa, S, yerr = S_std, color = 'green', linestyle='--', label = label1, lw = lw, ms = ms)
        ax.errorbar(aa, H2, yerr = H2_std, color = 'orange', linestyle='-', label = label2, lw = lw, ms = ms)
        ax.errorbar(aa, S2, yerr = S2_std, color = 'green', linestyle='-', label = label2, lw = lw, ms = ms)
        
    ax.set_xlim(0,len(aa)+1)
    
	#Another axis, just for having the sequence in ONE letter amino acid code
    ax2 = ax.twiny()
    ax2.plot(aa,H,c='none',lw=0.0)
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(aa)
    ax2.set_xticklabels(aa_name, Fontsize = SMALL_SIZE)
    # ax2.spines['bottom'].set_position(ax.spines['bottom'].get_position())
    ax2.tick_params(axis = u'both', which = u'both', length = 0, color = 'w') #remove tick lines
    # plt.setp(ax.get_xticklabels(), fontsize=SMALLEST_SIZE)
    
	#Make certain letters (the p-sites) in red so they stand out
    ax3 = ax.twiny()
    ax3.plot(aa,H,c='none',lw=0.0)
    ax3.set_xlim(ax.get_xlim())
    ax3.set_xticks(np.array(aa)[psites])
    ax3.set_xticklabels(np.array(aa_name)[psites], Fontsize = SMALL_SIZE, color = 'r')
    # ax2.spines['bottom'].set_position(ax.spines['bottom'].get_position())
    ax3.tick_params(axis = u'both', which = u'both', length = 0, color = 'w') #remove tick lines

    ax.set_xlabel("Res. Num.")
    ax.set_ylabel(r"% Helix")
    ax.set_ylim(0, 1)
    #ax.title(title)
    ax.legend(loc='upper right', bbox_to_anchor=(1.75, 1.0), fontsize = 8)
    # ax.set_title(title)
    # plt.savefig('./Figures/SSP-TraDES-BME.png', dpi=300, transparent=False, bbox_inches='tight')
	
	#save in a folder in the current directory called figures ('./Figures/') with the title you want. Can change .png to .pdf if you like.
    plt.savefig(path + 'dssp_' + ss + "_" + title + '.png', dpi=900, transparent=False, bbox_inches='tight')
    plt.savefig(path + 'dssp_' + ss + "_" + title + '.pdf', dpi=900, transparent=False, bbox_inches='tight')
    
    return aa

def trianalysis(top, traj, idx_start_arr, idx_end_arr, sample, num_protein, start_frame, end_frame, dt, residues):
    '''
    

    Parameters
    ----------
    top : top file
        input with mdtraj
    traj : trajectory
        input with mdtraj.
    idx_start_arr : array 
       array of atom start indeces
    idx_end_arr : array
        array of atom end indeces
    sample : string
        name of sample
    num_protein : int
        number of protein in system
    start_frame : int
        start frame for analysis
    end_frame : int
        end frame for analysis
    dt :int
        time step for analysis
    residues : list
        list of residues for analysis


    Returns
    -------
    H_arr : array of float
        Avg Helical propensity over residues per frame
    C_arr : array of float
        Avg Coil propensity over residues per frame
    S_arr : array of float
        Avg Strand propensity over residues per frame
    Ree_arr : array of float
        End-to-end distance
    Rg_arr : array of float
        Radius of gyration
    sample_arr : array of strings
        Sample name

    '''
    
    for i in range(0, num_protein):
        print("Analyzing protein %s" %(i+1))
            
        prot = top.select("index %s to %s and residue 0 to 39" %(idx_start_arr[i], idx_end_arr[i]))
        traj_prot = traj.atom_slice(prot)
        
        #Calculate avg EE and Rg for each frame
        nterm = top.select("index %s to %s and name CA and resid 0" %(idx_start_arr[i], idx_end_arr[i]))
        cterm = top.select("index %s to %s and name CA and resid 36" %(idx_start_arr[i], idx_end_arr[i]))
        atom_arr = np.array([[nterm[0], cterm[0]]])
        dist = md.compute_distances(traj_prot, atom_arr)
        
        #Calculate Rg for each frame
        rg = md.compute_rg(traj_prot)

        C_arr = []
        S_arr = []
        H_arr = []
        
        Ree_arr = []
        Rg_arr = []
        sample_arr = []
        
        t = 0
        for frame in range(start_frame, end_frame, dt):
            
            coil=np.zeros(traj_prot.topology.n_residues)
            helix=np.zeros(traj_prot.topology.n_residues)
            strand=np.zeros(traj_prot.topology.n_residues)
            
            print("Frame %s" %(frame))
        
            dssp=md.compute_dssp(traj_prot[frame])[0]
            for res in range(len(dssp)):
                if 'C' in dssp[res]:
                    coil[res]=1
                elif 'H' in dssp[res]:
                    helix[res]=1
                else:
                    strand[res]=1
                  
            #Update time 
            t+=1
                
            helix_avg = np.mean(helix[residues])
            coil_avg = np.mean(coil[residues])
            strand_avg = np.mean(strand[residues])
                
            H_arr.append(helix_avg)
            C_arr.append(coil_avg)
            S_arr.append(strand_avg)
            
            Ree_arr.append(dist[frame][0])
            Rg_arr.append(rg[frame])
            sample_arr.append(sample)
            
            ### Marked for removal ###
            #traj_slice = traj.atom_slice(prot)[frame:frame+1]
            #C = coil, #S = strand, H = helix
            
            # C1_arr.append(C1)
            # S1_arr.append(S1)
            
    return H_arr, C_arr, S_arr, Ree_arr, Rg_arr, sample_arr
