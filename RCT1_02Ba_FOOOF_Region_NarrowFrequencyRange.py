#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:42:14 2025

@author: riannehaartsen
"""

###################################################################################################
# FOOOF in RCT1 arbaclofen study: for narrow frequency range per region per condition
# --------------------------------------
#

###################################################################################################
# General imports
import numpy as np
import pandas as pd

# Import the FOOOF 
from fooof import FOOOF
# Import other functions
from matplotlib import pyplot as plt
import glob
#import csv


# create a list of file names for the loop
path = 'xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region'
powspec_names = glob.glob(path + '/*.csv')

print(len(powspec_names))

a, b = '.csv', ''
for i, v in enumerate(powspec_names):
    if a in v:
        powspec_names[i] = v.replace(a, b)
        

a, b = 'xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/', ''
for i, v in enumerate(powspec_names):
    if a in v:
        powspec_names[i] = v.replace(a, b)
        
        
print(powspec_names)

with open('xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/powspec_names.txt', 'w') as f:
    for item in powspec_names:
        f.write("%s\n" % item)

#%% open the file with names

fnamefile = open("xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/powspec_names.txt", "r")

# frequencies for power spectra data
freqs = np.arange(0, 120.5, 0.5).tolist()

#%%

all_res = {}

for aline in fnamefile.readlines():
    alldata = []
    data = pd.read_csv('xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/' + aline.rstrip() + '.csv')
    
    alldata.append(data) #reappend the new columns to the old data file

    pid = aline.rstrip()
    
    freq=(np.array(freqs))
    
    spectrum_x=(np.array(data))
    spectrum = np.reshape(spectrum_x,[241, ])

   
    fm=FOOOF(peak_width_limits=[1, 8], max_n_peaks=4, peak_threshold=0.1, aperiodic_mode='fixed') 
    
    freq_range=[3,28] 
    
    fm.fit(freq, spectrum, freq_range)
    
    ap_params, peak_params, r_squared, fit_error, gauss_params = fm.get_results()
    
    fm.report(freq, spectrum, freq_range)
    
    fm.plot()
    
    plt.savefig('xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/Narrow_FBoI/Narrow_' + aline.rstrip()  + '.png')

    plt.clf()
    
    plt.close()

    pp = peak_params[0:3,0] 
    
    num_peaks=len(pp)
    
    if len(pp) < 1:
        p1=0
        p2=0
        p3=0
    
    elif len(pp) == 1:
        p1= pp[0]
        p2=0
        p3=0
        
    elif len(pp) == 2:
        p1, p2 = pp[0], pp[1]
        p3=0
        
    elif len(pp) == 3:
        p1, p2, p3 = pp[0], pp[1], pp[2]  
        
    amp = peak_params[0:3,1]
    
    if len(amp) < 1:
        a1=0
        a2=0
        a3=0
    
    elif len(amp) == 1:
        a1= amp[0]
        a2=0
        a3=0
        
    elif len(amp) == 2:
        a1, a2 = amp[0], amp[1]
        a3=0
        
    elif len(amp) == 3:
        a1, a2, a3 = amp[0], amp[1], amp[2]  
    
    all_res[pid] = [r_squared,  fit_error, ap_params[0], ap_params[1], num_peaks, p1, p2, p3, a1, a2, a3]

    export = pd.DataFrame.from_dict(all_res, orient='index',)
    export.rename(index=str, columns={0: "r squared", 1: "fit error", 2:"intercept", 3:"slope", 4:'number of peaks', 
                                      5: 'peak 1 freq', 6:'peak 2 freq', 7:'peak 3 freq', 8:'peak 1 amplitude',
                                      9: 'peak 2 amplitude', 10: 'peak 3 amplitude'}, inplace=True)
    
    export.to_csv('xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/Narrow_FBoI/RCT1_Region_Narrow_data.csv')    
    