#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:25:55 2023

@author: tombullock

PHYSIO_TREV_ANALYSIS

"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# set directories
parentDir = '/Users/tombullock/Documents/Psychology/TREV_DEV'
dataDir = os.path.join(parentDir,'Data')
scriptDir = os.path.join(parentDir, 'TREV_Pipeline')
outDir = dataDir # for now

# set subject(s) - just do a single one for now
subjects = [10,11,12,13,15,16]

# LOOP THROUGH SUBJECTS
d=[]
for sjNum in subjects:
    
    # load TREV data
    filename = f'sj{sjNum:02d}_RADARS.csv'
    filepath = os.path.join(dataDir,filename)
    df_cardiac = pd.read_csv(filepath)
    df_cardiac.rename(columns={df_cardiac.columns[0]:'beat'},inplace=True) # rename beat index column
    
    # load event (trigger) data
    filename = f'sj{sjNum:02d}_RADARS_events.csv'
    filepath = os.path.join(dataDir,filename)
    df_event = pd.read_csv(filepath)
    
    # load bp data
    filename = f'sj{sjNum:02d}_RADARS_bp.csv'
    filepath = os.path.join(dataDir,filename)
    df_bp = pd.read_csv(filepath)
        
    # loop through all blocks, parse relevant data, do calculations, export etc.
    events_task = ['11', '12', '13', '21', '22', '23', '31','32','33','41','42','43', '51','52','53']
    events_rest = ['61','62','63','64','65','81','82','83','84','85','58']
    events_all = events_task #+ events_rest [events rest is causing problems so ignore for now]
    contractility_task = []
    contractility_rest = []
    for event in events_all:
        
        event = int(event)
        
        # get event timestamp for start block    
        timestamp_start = df_event[df_event.code==event]['time'].item()
        
        # find start idx of closest heartbeat in cardiac data file
        idx_start = np.argmin(abs(df_cardiac.time - timestamp_start))
        
        # get event timestamp for end block
        if event<55:
            block_duration = 90
            timestamp_end = timestamp_start + block_duration # task blocks are 90 s dur
        else:
            block_duration = 210
            timestamp_end = timestamp_start + block_duration # rest blocks are 210 s dur
            
        # find end idx of closest heartbeat in cardiac data file
        idx_end = np.argmin(abs(df_cardiac.time - timestamp_end))
        
        # grab peak times and contractility from data chunk and stick into dictionary
        
        # adjust peak times so zero is start of block
        peak_times = df_cardiac.iloc[idx_start:idx_end]['time'].values
        peak_times = peak_times-peak_times[0]
        
        contractility = df_cardiac.iloc[idx_start:idx_end]['contractility'].values
        
        n_beats = len(df_cardiac.iloc[idx_start:idx_end]['time'].values)
        
        x = np.linspace(1,block_duration,block_duration)
        
        # interpolate contractility so that each task block is 90s and each rest block is 210s
        contractility_interp = np.interp(x,peak_times,contractility)
        
        ## find start/end indices for closest BP measure in BP data file then grab bp data
        idx_start_bp = np.argmin(abs(df_bp.time - timestamp_start))
        idx_end_bp = idx_start_bp+block_duration
        bp = df_bp.iloc[idx_start_bp:idx_end_bp]['bp'].values

        # separate by event code into meaningful conditions/sessions
        if event in np.array([11,21,31,41,51]):
            condition='Pre'
        elif event in np.array([12,22,32,42,52]):
            condition='CPT'
        elif event in np.array([13,23,33,43,53]):
            condition='Post'
        
        # get session code from event
        session = 'S' + str(event)[0]
            
        # create a dict
        d.append({
            'sjNum':sjNum,
            'condition':condition,
            'session':session,
            'beats_count': n_beats,
            'hr_mean':round((n_beats/len(x))*60,1), # check this
            'contractility': contractility_interp,
            'contractility_mean':contractility_interp.mean(),
            'timepoints': x,
            'bp_mean':bp.mean(),
            'bp_median':np.median(bp),
            'bp':bp
            })
        
        del timestamp_start, timestamp_end, idx_start, idx_end, block_duration
        
# create a dataframe for (d)
df = pd.DataFrame(d)

# export dataframe to file
filename = 'PHYSIO_COMPILED.csv'
filepath = os.path.join(outDir,filename)
df.to_csv(filepath)







# create a quick boxplot for HR averaged across block
# sns.boxplot(
#     data=df,
#     x='condition',
#     y='heartrate',
#     hue='session'
#     )

# quick boxplot for hr_mean
sns.boxplot(
    data=df,
    x='session',
    y='hr_mean',
    hue='condition'
    )
    
# # quick boxplot for bp_mean
# sns.boxplot(
#     data=df,
#     x='session',
#     y='bp_mean',
#     hue='condition'
#     )

# # quick boxplot for contractility_mean
# sns.boxplot(
#     data=df,
#     x='session',
#     y='contractility_mean',
#     hue='condition'
#     )


# # plot contractility averaged over sessions for each cond
# cont_pre = (contractility_task[0]+contractility_task[3]+contractility_task[6]+contractility_task[9]+contractility_task[12])/5
# cont_cpt = (contractility_task[1]+contractility_task[4]+contractility_task[7]+contractility_task[10]+contractility_task[13])/5
# cont_post = (contractility_task[2]+contractility_task[5]+contractility_task[8]+contractility_task[11]+contractility_task[14])/5

# data = cont_pre #d[0]['contractility_interp']
# timepoints = timepoints_task #d[0]['timepoints']
# sns.lineplot(x=timepoints,y=data)
    
# data = cont_cpt # d[1]['contractility_interp']
# timepoints = timepoints_task #d[1]['timepoints']
# sns.lineplot(x=timepoints,y=data)
        
# data = cont_post #d[2]['contractility_interp']
# timepoints = timepoints_task #d[2]['timepoints']
# sns.lineplot(x=timepoints,y=data)


# plt.legend(['pre','cpt','post'])

# plt.title("Contractility")
    
# plt.xlabel('Time(s)')    
    
   
    
    
