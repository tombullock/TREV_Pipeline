#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:26:11 2023

@author: neil (tom edits to add events)
"""

import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilename
import pandas as pd
from scipy.signal import find_peaks as findpeaks
from scipy.signal import firwin, lfilter, filtfilt, butter
from scipy.stats import zscore
import bioread
import os
import matplotlib.pyplot as plt




# CELL 1
from miniMEAPsys import loadem
cont_dict,file_path=loadem()
print('data loaded!')


# CELL 2
from miniMEAPsys import compute_peaks,thresh_ind
%matplotlib tk
plt.close('all')
i=thresh_ind(cont_dict['hz'])
plt.plot(cont_dict['t'][i],cont_dict['s'][i])
peak_times,peak_vals,peak_inds=compute_peaks(cont_dict)

plt.close('all')



## CELL 5: OUTPUT DATA TO CSV ###
# this cell will model respiration and heart-rate out of the beatwise contractility estimates and output to a csv
import pandas as pd
from miniMEAPsys import regress_out,outsave_box
plt.close('all')
cont_dict['cont_peak_times']=peak_times
cont_dict['raw_contractility_peaks']=peak_vals
cont_dict=regress_out(cont_dict)

out_f = outsave_box(file_path[:-4]+'.csv')

out_dict={}
out_dict['time']=cont_dict['cont_peak_times']
out_dict['contractility']=cont_dict['resid_contractility']
pd.DataFrame(out_dict).to_csv(out_f)

# TOM ADDED EVENT CODES OUTPUT TO SEPARATE FILE
if 'events_s' in cont_dict: # dict will only contain "events_s" key if event channel present and selected in prepro
    out_f_event = out_f[:-4] + '_events.csv' #outsave_box(file_path[:-4]+'_events.csv')
    
    events = cont_dict['events_s']
    times = cont_dict['events_t']
    
    code_mat = []
    time_mat = []
    
    # cut out redundant events from oversampling
    for i in np.arange(np.shape(events)[0]): 
        if events[i] != 0 and events[i-1] != events[i]:       
            code_mat.append(int(events[i]))
            time_mat.append(round(times[i],3))
    
    # create df with event codes and timestamps and write to csv
    pd.DataFrame({'time':time_mat,'code':code_mat}).to_csv(out_f_event)

