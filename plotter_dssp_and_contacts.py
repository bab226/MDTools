#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 12:58:36 2023

Before using:
Use dssp_rg_and_ree_mdtraj.py and contact_prob_PDF.py to analyze data.

Description: This script will unpickle data and plot results.

@author: bab226
"""
import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import random 
from matplotlib.ticker import FormatStrFormatter

# COLOR

def get_cmap(n, name='viridis'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

color = sns.color_palette('deep', 5, desat=1)
color2 = sns.color_palette('coolwarm', 20, desat=1)



# TEXT

matplotlib.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['mathtext.default'] = "regular"
matplotlib.rcParams['axes.linewidth'] = 1 # set the value globally

SMALLEST_SIZE = 2.0 #was 1.5
SMALL_SIZE = 12/2.2
MEDIUM_SIZE = 18/2.25
BIGGER_SIZE = 20/2.25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



# BODY
title = "IAPP_ADM"
path = '/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/'

#Load biophysical information:
woADM = open(path + "dssp_wo_adm_100-200ns.pkl", "rb")    # rb is reading binary data
wADM = open(path + "dssp_adm_100-200ns.pkl", 'rb')

woADM = pickle.load(woADM)
wADM = pickle.load(wADM)
aa = wADM['resname']

#Load probability density distribution for drug binding:
PDF = np.load('/Users/bab226/Documents/yale_research/iapp/md_sims/remd_iapp/iapp_adm116/new_files/analysis/contacts/contact-density.npz')
residues = PDF['a']
prob = PDF['b']



#Plot data with SS from with and without drug with binding profile of inhibitor.

width = 3.42
fig, ax1 = plt.subplots(1, 1, figsize=(width,width/1.2))
plt.subplots_adjust(left=0.2, right=0.85, top=0.83, bottom=0.18)

ax1.axvspan(5, 15, alpha=0.1, color = color2[0], lw=None)
ax1.axvspan(23, 27, alpha=0.1, color = color2[0], lw=None)

ax1.plot(woADM['helix'], linestyle='--', color = color[2], label = "IAPP")
ax1.set_ylim([0,1])
ax1.set_xlim(0,len(aa)+1)

ax1.plot(wADM['helix'], linestyle='-', color = color[1],  label = "IAPP + SM")
ax1.set_ylim([0,1])
ax1.set_xlabel("Residue #")
ax1.set_ylabel("Fraction Helicity")
ax1.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fontsize = 8)
#Another axis, just for having the sequence in ONE letter amino acid code.


ax2 = ax1.twiny()
ax2.plot(aa, wADM['helix'],c='none',lw=0.0)   #Plot invisible data set
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(residues)   #Add ticks every residue
ax2.set_xticklabels(aa, Fontsize = SMALL_SIZE)  #Label each tick with AA One Letter Code
ax2.tick_params(axis = u'both', which = u'both', length = 0, color = 'w') #remove tick lines



plt.savefig(path + 'dssp_' + title + '.png', dpi=900, transparent=False, bbox_inches='tight')

width = 3.42
fig, ax1 = plt.subplots(1, 1, figsize=(width,width/1.2))
plt.subplots_adjust(left=0.2, right=0.85, top=0.83, bottom=0.18)

ax1.axvspan(5, 15, alpha=0.1, color = color2[0], lw=None)
ax1.axvspan(23, 27, alpha=0.1, color = color2[0], lw=None)

ax1.plot(residues, prob, linestyle='-', color = color[0], label = "IAPP")
ax1.set_xlim(ax1.get_xlim())
ax1.set_ylim([0, 0.1])
ax1.set_xlim(0,len(aa)+1)
ax1.set_xlabel("Residue #")
ax1.set_ylabel("PDF")

#Another axis, just for having the sequence in ONE letter amino acid code.

ax2 = ax1.twiny()
ax2.plot(aa, wADM['helix'],c='none',lw=0.0)   #Plot invisible data set
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(residues)   #Add ticks every residue
ax2.set_xticklabels(aa, Fontsize = SMALL_SIZE)  #Label each tick with AA One Letter Code
ax2.tick_params(axis = u'both', which = u'both', length = 0, color = 'w') #remove tick lines
plt.savefig(path + 'contacts_PDF_' + title + '.png', dpi=900, transparent=False, bbox_inches='tight')