#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 15:23:27 2023

@author: tombullock

PHYSIO_Parse_TREV_Prepro_Data

Purpose: 

Load data outputted from TREV prepro pipeline (requires sjXX_RADARS.csv and sjXX_RADARS_events.csv).
Parse the beat timing and contractility data into chunks based on events and save into a dict (d)
ToDo: Compute cardaic measures

"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os
import janitor


# set dirs 
rDir = '/Users/tombullock/Documents/Psychology/RADARS'
sourceDirData = 'Data_TREV'
#sourceDirImages = 'Plots'
#destDir = 'Plots_Py'
destDirData = 'Data_Compiled_Py'

# subject for processing
sjNum = 10

# load preprocessed TREV data + events file
trev_filepath = os.path.join(rDir,sourceDirData,'sj'+str(sjNum)+'_RADARS.csv')
events_filepath = os.path.join(rDir,sourceDirData,'sj'+str(sjNum)+'_RADARS_events.csv')

df_trev = pd.read_csv(trev_filepath)
df_events = pd.read_csv(events_filepath)


## Parse files according to events

# create lists of codes (here's where you'd add subject specific exceptions if needed e.g. subject missing a rest block)
task_block_codes = [11,12,13,21,22,23,31,32,33,41,42,43,51,52,53]
rest_block_codes = [58,61,62,63,64,65,81,82,83,84,85]

# loop through the different trial types
d=[]
for trial_type in range(2): 
    
    if trial_type==0: # task blocks
        block_codes = task_block_codes
        block_duration = 90 # secs
    elif trial_type==1: # rest blocks
        block_codes = rest_block_codes
        block_duration = 210 # secs
        
    for this_code in block_codes:
        
        # find start and end of block based on event code start + 90s duration
        start_time = df_events[df_events['code']==this_code]['time'].item()
        start_idx_trev_file = abs(df_trev['time']-start_time).idxmin()
        #start_time = df_trev['time'].loc[start_idx_trev_file]
        end_time = start_time + block_duration
        end_idx_trev_file = abs(df_trev['time']-end_time).idxmin()
        
        # get arrays of beat times and contractility
        beat_times = df_trev['time'].loc[start_idx_trev_file:end_idx_trev_file].values
        beat_times = beat_times - beat_times[0] # reset so first beat is at zero
        contractility = df_trev['contractility'].loc[start_idx_trev_file:end_idx_trev_file].values
        
        # compute cardiac measures here! [Caitlin?]
        
        # append measures to new dict/df [probably want to add the various cardiac measures to this dict too]
        d.append({
            'Code':this_code,
            'Peak_Times':beat_times,
            'Contractility':contractility          
            })
        
        # note if needed you can access these in the dict with 
        # e.g. d[0]['Peak_Times']
        
        del(start_time,start_idx_trev_file,end_time,end_idx_trev_file,beat_times,contractility)
    
    
#df = pd.DataFrame(d)
    
    
    
    
        
    


