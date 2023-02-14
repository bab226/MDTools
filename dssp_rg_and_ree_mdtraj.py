#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:30:04 2022

Mdanalysis script to analyze Rg vs. Rh in high throughput manner.

Adapted from Dr. Gregory 'Beggory-Sherlock' Gomes

MDTraj Citation:
@article{McGibbon2015MDTraj,
    author = {McGibbon, Robert T. and Beauchamp, Kyle A. and Harrigan, Matthew P. and Klein, Christoph and Swails, Jason M. and Hern{\'a}ndez, Carlos X.  and Schwantes, Christian R. and Wang, Lee-Ping and Lane, Thomas J. and Pande, Vijay S.},
    title = {MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories},
    journal = {Biophysical Journal},
    volume = {109},
    number = {8},
    pages = {1528 -- 1532},
    year = {2015},
    doi = {10.1016/j.bpj.2015.08.015}
}

@author: bab226
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import math
import pandas as pd
import mdtraj as md
from scipy import stats

from random import seed
from random import randint
from random import random
from random import uniform

from contact_map import ContactFrequency, ContactDifference, OverrideTopologyContactDifference

import sys
sys.path.append("/Users/bab226/Documents/yale_research/iapp/md_sims/house_scripts/mdanalysis/")
import mdtraj_helper as mdh

# Text size for plots.
SMALLEST_SIZE = 2.0 #was 1.5
SMALL_SIZE = 12/2.25
MEDIUM_SIZE = 18/2.25
BIGGER_SIZE = 20/2.25

# Font for labels of plots.
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Path to trajectory and topology files.
path = '/Users/bab226/Documents/yale_research/iapp/md_sims/iapp_adm116/new_files/analysis/'
aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


#Define protein 1 atom selections
num_protein = 1
num_atoms = 537
idx_start = 0
idx_end = idx_start + (num_atoms-1)
text1 = "a %d-%d" %(idx_start, idx_end)
idx_start_arr = []
idx_end_arr = []
idx_start_arr.append(idx_start)
idx_end_arr.append(idx_end)
#print(text1)
for i in range(1, num_protein+1):
    idx = i-1
    
    idx_start = idx_end + 1
    idx_end = idx_start + (num_atoms - 1)
    idx_start_arr.append(idx_start)
    idx_end_arr.append(idx_end)
    text1 = "a %d-%d" %(idx_start, idx_end)
    #print(text1)
    
#Define protein 2 atom selections (if using 2 different topology files)
#idx2_start = 2604
idx2_start = 0
idx2_end = idx2_start + (num_atoms-1)
text2 = "a %d-%d" %(idx2_start, idx2_end)
idx2_start_arr = []
idx2_end_arr = []
idx2_start_arr.append(idx2_start)
idx2_end_arr.append(idx2_end)
#print(text2)
for i in range(1, num_protein+1):
    idx2 = i-1
    
    idx2_start = idx2_end + 1
    idx2_end = idx2_start + (num_atoms - 1)
    idx2_start_arr.append(idx2_start)
    idx2_end_arr.append(idx2_end)
    text2 = "a %d-%d" %(idx2_start, idx2_end)
    #print(text2)

    
# Load Traj and topology files for MDAnalysis.
path = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/"

top1 = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_monomer/analysis/final.gro"
file1 = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_monomer/analysis/final_100-200ns.xtc"

top2 = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/final.gro"
file2 = "/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/final_100-200ns.xtc"


# Load Traj and topology files for MDTraj

mdtop1 = md.load(top1).topology
mdtraj1 = md.load(file1,top=top1,stride=1,)

mdtop2 = md.load(top2).topology
mdtraj2 = md.load(file2,top=top2,stride=1,)


#%%
''' Used to average/bin frames, but they cannot be correlated since adjacent 
structures are unrelated to eachother.

Now calculates Rg and Ree for each frame. Can delete commented code. '''
#for i in range(0, len(idx_start_arr)):
'''
for i in range(0, num_protein):
    print("Analyzing protein %s" %(i+1))
    bb = u1.select_atoms('protein and backbone and index %s-%s' %(idx_start_arr[i], idx_end_arr[i])) # a selection (AtomGroup)

    Ree=[]
    Rg=[]
    
    Rg_avg = []
    Rg_std = []
    
    Ree_avg = []
    Ree_std = []
    
    A=[]
    time=[]
    bins = 10
    
    t1 = 0
    t2 = 5000
    
    width = int(len(u1.trajectory[t1:t2])/bins)-1  # subtract one since starts from zero
    
    start = 0
    end = start + width
    
    bin1 = []
    bin2 = []
    
    print(start)
    print(end)
    bin1.append(start)
    bin2.append(end)
    for i in range(0, bins-1):
        start = end + 1
        end = start + width
        print(start)
        print(end)
        bin1.append(start)
        bin2.append(end)
        
        for ts in u1.trajectory[start:end]:     # iterate through all frames
            
            nterm = bb.select_atoms('protein and name CA')[0]
            cterm = bb.select_atoms('protein and name CA')[-1]
            
            r = cterm.position - nterm.position # end-to-end vector from atom positions
            # print(np.shape(bb[:].position))
            
            Ree.append(np.linalg.norm(r))  # end-to-end distance
            Rg.append( bb.radius_of_gyration())  # method of AtomGroup
            A.append(bb.asphericity())
            time.append(ts.frame)
    
        start = start + i*width
        
        Ree_avg.append(np.mean(Ree))
        Ree_std.append(np.std(Ree))
        
        Rg_avg.append(np.mean(Rg))
        Rg_std.append(np.std(Rg))
    
    width = 3.42
    
    fig, ax = plt.subplots(1, 1, figsize=(width,width/1.2))
    plt.xlabel('Frame')
    plt.ylabel(r'$R_{ee}$')
    plt.errorbar(bin1[:-1], Ree_avg, yerr=Ree_std, elinewidth=2)
    #plt.ylim(40,90)
    
    fig2, ax = plt.subplots(1, 1, figsize=(width,width/1.2))
    plt.xlabel('Frame')
    plt.ylabel(r'$R_{g}$')
    plt.errorbar(bin1[:-1], Rg_avg, yerr=Rg_std, elinewidth=2)
    #plt.ylim(6,30)
    
    
#%%
# g1 = sns.jointplot(x=Rg, y=Ree,kind="hex", color="g")
# # or set labels via the axes objects
# g1.ax_joint.set_xlabel(r'$R_{g}$', fontweight='bold')
# g1.ax_joint.set_ylabel(r'$R_{ee}$', fontweight='bold')

# g = sns.jointplot(x=Rg, y=A,kind="hex", color="g",)
# # or set labels via the axes objects
# g.ax_joint.set_xlabel(r'$R_{g}$', fontweight='bold')
# g.ax_joint.set_ylabel(r'$A$', fontweight='bold')


# g = sns.jointplot(x=Ree, y=A,kind="hex", color="g")
# # or set labels via the axes objects
# g.ax_joint.set_xlabel(r'$R_{ee}$', fontweight='bold')
# g.ax_joint.set_ylabel(r'$A$', fontweight='bold')

'''
#bins = 1

t1 = 0
t2 = len(mdtraj1)

dt = 1   #Timestep for analysis

#width = int(np.round((len(mdtraj1))/bins))

start = 0

sample = ['wo-ADM', 'w-ADM']
sample_arr = []

residues = list(range(3,27))


#end = start + width
# bin1 = []
# bin2 = []

# bin1.append(start)
# bin2.append(end)

# if bins > 1:
#     for i in range(0, bins-1):
#         start = end + 1
#         end = start + width-1
#         print(start)
#         print(end)
#         bin1.append(start)
#         bin2.append(end)
#     else: 
#         bin1.append(0)
#         bin2.append(len(mdtraj1))
 
     #Calculate Ree and Rg with mdtraj.
     #AVERAGED FRAMES are uncorrelated!
     #Calculate avg EE and Rg for each frame and compare to DSSP
     
H1_arr, C1_arr, S1_arr, Ree_arr, Rg_arr, sample_arr = mdh.trianalysis(mdtop1, mdtraj1, 
        idx_start_arr, idx_end_arr, sample[0], num_protein, t1, t2, dt, residues)

H2_arr, C2_arr, S2_arr, Ree_arr2, Rg_arr2, sample_arr2 = mdh.trianalysis(mdtop2, mdtraj2, 
        idx_start_arr, idx_end_arr, sample[1], num_protein, t1, t2, dt, residues)



#%%
# Create Pandas dataframe to hold relevent data.
data1 = {'Rg': Rg_arr,
        'Ree': Ree_arr,
     #   'RgSD': Ree_std,
     #   'ReeSD': Ree_std,
        'helicity' : H1_arr,
        'coil' : C1_arr,
        'strand' : S1_arr,
     #   'helicitySD' : mean_helix_std,
        'Sample': sample_arr}

df1 = pd.DataFrame(data1)

data2 = {'Rg': Rg_arr2,
        'Ree': Ree_arr2,
     #   'RgSD': Ree_std,
     #   'ReeSD': Ree_std,
        'helicity' : H2_arr,
        'coil' : C2_arr,
        'strand' : S2_arr,
     #   'helicitySD' : mean_helix_std,
        'Sample': sample_arr2}

df2 = pd.DataFrame(data2)

frames = [df1, df2]
all_data = pd.concat(frames)

df = pd.DataFrame(all_data)
#%% Save data in df
df1.to_pickle(path + "data_wo_adm_100-200ns.pkl")
df2.to_pickle(path + "data_adm_100-200ns.pkl")
df.to_pickle(path + "data_100-200ns.pkl")

#%%
# 2D Plots for correlation detection.
#################################
wdth = 3.42*1.2
#Rg vs. Ree
df_filter = df

g1 = sns.jointplot(data=df_filter, x=df_filter['Rg'], y=df_filter['Ree'],
                    height=wdth, ratio=3,  
                    #xlim=[1, 1.25], ylim=[-0.01, 0.04], 
                    hue="Sample", kind='hist')
                    #joint_kws=dict(bw_adjust = 0.5, levels=10))
# g1 = sns.jointplot(data=df, x=df['Rg']/10, y=df['Ree']/10, color="g",  
#                   height=width, ratio=3, s=90, hue="Sample")
ax = plt.gca()
sns.move_legend(ax, "upper left")

# or set labels via the axes objects
g1.ax_joint.set_xlabel(r'$R_{g}$ (nm)', fontweight='normal', fontsize=22)
g1.ax_joint.set_ylabel(r'$R_{ee}$ (nm)', fontweight='normal', fontsize=22)

plt.savefig(path + "ree_vs_rg.png", dpi=600, transparent=False, bbox_inches='tight')
plt.savefig(path + "ree_vs_rg.pdf", dpi=600, transparent=False, bbox_inches='tight')

#################################

#Rg vs. helicity
g2 = sns.jointplot(data=df_filter, x=df_filter['Rg'], y=df_filter['helicity'],
                    height=wdth, ratio=3, 
                    #xlim=[1, 1.25], ylim=[-0.01, 0.04], 
                    hue="Sample", kind='hist')
                    #joint_kws=dict(bw_adjust = 0.5, levels=10))
ax = plt.gca()
sns.move_legend(ax, "upper left")

# or set labels via the axes objects
g2.ax_joint.set_xlabel(r'$R_{g}$ (nm)', fontweight='normal', fontsize=22)
g2.ax_joint.set_ylabel('mean helicity', fontweight='normal', fontsize=22)

plt.savefig(path + "rg_vs_helicity.png", dpi=600, transparent=False, bbox_inches='tight')
plt.savefig(path + "rg_vs_helicity.pdf", dpi=600, transparent=False, bbox_inches='tight')

#################################

# Helicity vs. Ree
g3 = sns.jointplot(data=df_filter, x=df_filter['helicity'], y=df_filter['Ree'],
                    height=wdth, ratio=3, 
                    #xlim=[-0.01, 0.04],ylim=[1.8,2.5],
                    hue="Sample", kind='hist')

ax = plt.gca()
sns.move_legend(ax, "upper left")

# or set labels via the axes objects
g3.ax_joint.set_xlabel('mean helicity', fontweight='normal', fontsize=22)
g3.ax_joint.set_ylabel(r'$R_{ee}$ (nm)', fontweight='normal', fontsize=22)

plt.savefig(path + "helicity_vs_ree.png", dpi=600, transparent=False, bbox_inches='tight')
plt.savefig(path + "helicity_vs_ree.pdf", dpi=600, transparent=False, bbox_inches='tight')

# Helicity vs. Strand
g4 = sns.jointplot(data=df_filter, x=df_filter['helicity'], y=df_filter['strand'],
                    height=wdth, ratio=3, 
                    #xlim=[-0.01, 0.04],ylim=[1.8,2.5],
                    hue="Sample", kind='hist')

ax = plt.gca()
sns.move_legend(ax, "upper left")

# or set labels via the axes objects
g3.ax_joint.set_xlabel('mean helicity', fontweight='normal', fontsize=22)
g3.ax_joint.set_ylabel('mean strand', fontweight='normal', fontsize=22)

plt.savefig(path + "helicity_vs_strand.png", dpi=600, transparent=False, bbox_inches='tight')
plt.savefig(path + "helicity_vs_strand.pdf", dpi=600, transparent=False, bbox_inches='tight')
#%%
# Helicity vs. Disorder
g4 = sns.jointplot(data=df1, x=df1['helicity'], y=df1['coil'],
                    height=wdth, ratio=3, 
                    xlim=[-0.01, 0.8],
                    ylim=[0, 1],
                    #hue="Sample", 
                    joint_kws=dict(fit_reg=True),
                    kind='reg')
                    
sns.jointplot(data=df2, x=df2['helicity'], y=df2['coil'],
                    height=wdth, ratio=3, 
                    xlim=[-0.01, 0.8],
                    ylim=[0, 1],
                    #hue="Sample", 
                    joint_kws=dict(fit_reg=True),
                    kind='reg')

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(df1['helicity'], df1['coil'])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(df2['helicity'], df2['coil'])

print ("""
regression 1: y = %.3f x + %.3f, p-val = %.3f, std_err = %.3f
r1-squared: %.3f
regression 2: y = %.3f x + %.3f, p-val = %.3f, std_err = %.3f
r1-squared: %.3f"""
        %(slope1, intercept1, p_value1, std_err1, r_value1**2, slope2, intercept2, r_value2, std_err2, p_value2**2))
        
ax = plt.gca()
#sns.move_legend(ax, "upper left")

# or set labels via the axes objects
g3.ax_joint.set_xlabel('mean helicity', fontweight='normal', fontsize=22)
g3.ax_joint.set_ylabel('mean coil', fontweight='normal', fontsize=22)

plt.savefig(path + "helicity_vs_coil.png", dpi=600, transparent=False, bbox_inches='tight')
plt.savefig(path + "helicity_vs_coil.pdf", dpi=600, transparent=False, bbox_inches='tight')

#%% 
# Summary Text


#################################

m1_rg = np.mean(df1.Rg)
m1_ree = np.mean(df1.Ree)
m1_helicity = np.mean(df1.helicity)

m2_rg = np.mean(df2.Rg)
m2_ree = np.mean(df2.Ree)
m2_helicity = np.mean(df2.helicity)

m1_rg_std = np.std(df1.Rg)
m1_ree_std = np.std(df1.Ree)
m1_helicity_std = np.mean(df1.helicity)

m2_rg_std = np.std(df2.Rg)
m2_ree_std = np.std(df2.Ree)
m2_helicity_std = np.std(df2.helicity)

print("""
      Mean Rg for %s is %.3f +/- %.3f nm
      Mean Ree for %s is %.3f +/- %.3f nm
      Mean Helicity for %s is %.3f +/- %.3f
      Mean Rg for %s is %.3f +/- %.3f nm
      Mean Ree for %s is %.3f +/- %.3f nm
      Mean Helicity for %s is %.3f +/- %.3f
      """ %(sample[0], m1_rg, m1_rg_std, sample[0], m1_ree, m1_ree_std, sample[0], m1_helicity, m1_helicity_std,
      sample[1], m2_rg, m2_rg_std, sample[1], m2_ree, m2_ree_std, sample[1], m2_helicity, m2_helicity_std))
      
#%% Contact Maps for IAPP

prot1 = mdtop1.select("index %s to %s and residue 1 to 37" %(idx_start_arr[0], idx_end_arr[0]))
traj_prot1 = mdtraj1.atom_slice(prot1)
prot2 = mdtop2.select("index %s to %s and residue 1 to 37" %(idx2_start_arr[0], idx2_end_arr[0]))
traj_prot2 = mdtraj2.atom_slice(prot2)

#cutoff is max distance for contact, 0.4 = 4A
traj_contacts1 = ContactFrequency(traj_prot1, cutoff=0.4)
traj_contacts2 = ContactFrequency(traj_prot2, cutoff=0.4)
diff = traj_contacts2 - traj_contacts1

fig, ax = traj_contacts1.residue_contacts.plot()
plt.xlabel("Residue")
_ = plt.ylabel("Residue")
plt.savefig(path + "contact_map_iapp.pdf", dpi=600, transparent=False, bbox_inches='tight')
# setting the return of plt.ylabel to the dummy var _ is a way to swallow output

fig, ax = traj_contacts2.residue_contacts.plot()
plt.xlabel("Residue")
_ = plt.ylabel("Residue")
plt.savefig(path + "contact_map_iapp_adm.pdf", dpi=600, transparent=False, bbox_inches='tight')

fig, ax = diff.residue_contacts.plot()
plt.xlabel("Residue")
_ = plt.ylabel("Residue")
plt.savefig(path + "contact_map_diff.pdf", dpi=600, transparent=False, bbox_inches='tight')

#%% For ADM-IAPP interaction
'''A Prob. density function was generated elsewhere. Please disregard this code and delete in future verison'''
'''
import itertools

prot1 = mdtop2.select("index %s to %s and residue 1 to 37" %(idx_start_arr[0], idx_end_arr[0]))
traj_prot1 = mdtraj1.atom_slice(prot1)
prot2 = mdtop2.select("index 616 621 622 623 625 626 627 628 629 630 631 632 633 634 636 638")
traj_prot2 = mdtraj2.atom_slice(prot2)
r1 = list(range(1,38))
r2 = [38]
pairs = list(itertools.product(r1, r2))
#%%
#cutoff is max distance for contact, 0.4 = 4A

dist, pairs = md.compute_contacts(mdtraj2, pairs, ignore_nonprotein=False)

res = list(range(0,37))
sample_arr = []
dist_matrix = []
for j in res:
    dist_arr = []
    for i in range(0, len(dist)):
        sample = "res %s_%s" %(pairs[j-1][0], pairs[j-1][1])
        sample_arr.append(sample)
        
        print("Resid %s and %s" %(pairs[j-1][0], pairs[j-1][1]))
        fr_dist = dist[i][j-1]
        dist_arr.append(fr_dist)
    
    dist_matrix.append(dist_arr)
    
time = list(range(0,len(dist_matrix[13])))

wth = 2.3 *1.2 #originall *1
hgt = wth/1.25    

lw = 1
ms = 1 
colors = plt.cm.jet(np.linspace(0,1,len(dist_matrix)))
fig, ax = plt.subplots(figsize=(wth, hgt))
mean_arr = []
for i in range(0,37):
    mean = np.mean(dist_matrix[i])
    mean_arr.append(mean)
    ax.plot(time, dist_matrix[i], color=colors[i])
    ax.set_xlabel("Frame")
    ax.set_ylabel("Distance (nm)")
    ax.set_ylim(0, 10)
    '''
    