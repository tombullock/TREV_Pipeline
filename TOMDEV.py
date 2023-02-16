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

## 3,1,1,14


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




## CELL 3: INTERACTIVE PLOT - PEAKS IN ACCELERATION TIME SERIES###
# run this cell to open the interactive plot to visualise and amend peaks in the acceleration time series
# accleration peaks are more robust to noise and are better estimates of heartbeat times
# left click - add peak
# right click - remove peak
# hold m and left click - use moving average to estimate peak
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from miniMEAPsys import *

%matplotlib tk
plt.close('all')
fig,ax = plt.subplots()
plt.subplots_adjust(bottom=0.25)

w = ax.plot(cont_dict['t'],cont_dict['s'])
peak_plot, = ax.plot(peak_times,peak_vals ,'vm')


## RADARS: start custom code (TOM BULLOCK) - add event markers into plot so we can visualize start/end blocks

# isolate event codes and times from dict
events = cont_dict['events_s']
times = cont_dict['events_t']
    
code_mat = []
time_mat = []

# cut out redundant events from oversampling
for i in np.arange(np.shape(events)[0]): 
    if events[i] != 0 and events[i-1] != events[i]:       
        code_mat.append(str(int(events[i])))
        time_mat.append(round(times[i],3))
        
# find start of steady state blocks and then calculate end of blocks so we can insert end triggers
start_block_event_codes_ssvep = ['11', '12', '13', '21', '22', '23', '31','32','33','41','42','43', '51','52','53']
start_block_event_times_ssvep = []
end_block_event_times_ssvep = []
for i in np.arange(len(start_block_event_codes_ssvep)):
    event = start_block_event_codes_ssvep[i]
    this_time = time_mat[code_mat.index(event)]    
    start_block_event_times_ssvep.append(this_time)
    end_block_event_times_ssvep.append(this_time + 90)

end_block_event_codes_ssvep = ['911', '912', '913', '921', '922', '923', '931','932','933','941','942','943', '951','952','953']

# create list of start codes and times for resting state blocks
start_block_event_codes_rest = ['61','62','63','64','65','81','82','83','84','85','58']
start_block_event_times_rest = []
for i in np.arange(len(start_block_event_codes_rest)):    
    event = start_block_event_codes_rest[i]
    this_time = time_mat[code_mat.index(event)]
    start_block_event_times_rest.append(this_time)
    
# creae list of end codes and times for resting state blocks
end_block_event_codes_rest = ['71','72','73','74','75','91','92','93','94','95','59']
end_block_event_times_rest = []
for i in np.arange(len(end_block_event_codes_rest)):
    event = end_block_event_codes_rest[i]
    this_time = time_mat[code_mat.index(event)]
    end_block_event_times_rest.append(this_time)
    
# merge lists of codes and times    
start_block_codes = start_block_event_codes_ssvep + start_block_event_codes_rest
start_block_times = start_block_event_times_ssvep + start_block_event_times_rest
end_block_codes = end_block_event_codes_ssvep + end_block_event_codes_rest
end_block_times = end_block_event_times_ssvep + end_block_event_times_rest

# plot green lines for start block codes and red lines for end block codes
plt.vlines(x=start_block_times,ymin=-5,ymax=5,linestyles=':',colors = 'g',linewidth=4)
plt.vlines(x=end_block_times,ymin=-5,ymax=5,linestyles=':',colors = 'r',linewidth=4)


## RADARS: end custom code (TOM BULLOCK)


# Set the axis and slider position in the plot
bottom_pos = plt.axes([0.2, 0.1, 0.65, 0.03],facecolor = 'white')
scroll_slider = Slider(bottom_pos,'time', -1,cont_dict['t'][-1])

# Make a vertically oriented slider to control the amplitude
right_pos = plt.axes([.95, 0.25, 0.0225, 0.63])
y_zoom = Slider(right_pos,label="Y zoom",valmin=0,valmax=1,valinit=.5,orientation="vertical")
left_pos = plt.axes([.05, 0.25, 0.0225, 0.63])
x_zoom = Slider(left_pos,label="X zoom",valmin=0,valmax=1,valinit=.5,orientation="vertical")

ax.set_ylim(-np.median(peak_vals)*4*.5,np.median(peak_vals)*4*.5)

ax.set_xlim(-1,9)

orig_n = len(peak_times)

ax.set_title('Cell 3: remove / add acceleration peaks. Orig: '+str(orig_n)+', Updated: '+str(orig_n))





def update(val):
    pos = scroll_slider.val
    yzoom = y_zoom.val
    xzoom = x_zoom.val
    ax.set_xlim(pos, pos+20*xzoom)
    fig.canvas.draw_idle
    ecc=np.median(peak_vals)*4*yzoom
    ax.set_ylim(-ecc, ecc)
    fig.canvas.draw_idle
    
# update function called using on_changed() function
scroll_slider.on_changed(update)
y_zoom.on_changed(update)
x_zoom.on_changed(update)

# Display the plot
plt.show()

def onclick(event):
    global peak_times,peak_vals,fig,peak_plot,orig_n,ax
    ix, iy, ib, ik, iax = event.xdata, event.ydata, event.button, event.key, event.inaxes
    if (ix!=None) & (iax==ax):
        if (np.min(np.abs(ix-peak_times))<100) & (event.button!=1) & (ik==None):
            peak_times,peak_vals,peak_plot=remove_point(ix,peak_times,peak_vals,peak_plot)
            ax.set_title('Cell 3: remove / add acceleration peaks. Orig: '+str(orig_n)+', Updated: '+str(len(peak_vals)))
        elif (event.button==1) & (ik==None) & (ix>0) & (ix>0) & (iy<np.max(peak_vals)*1.2): #update dist from pos
            peak_times,peak_vals,peak_plot=add_point(ix,iy,peak_times,peak_vals,peak_plot)
            ax.set_title('Cell 3: remove / add acceleration peaks. Orig: '+str(orig_n)+', Updated: '+str(len(peak_vals)))
        elif (event.button==1)&(ik=='m'):
            peak_times,peak_vals,peak_plot=meap_dat(ix,iy,peak_times,peak_vals,peak_plot)
            ax.set_title('Cell 3: remove / add acceleration peaks. Orig: '+str(orig_n)+', Updated: '+str(len(peak_vals)))
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        if not plt.fignum_exists(1):
            canvas.mpl_disconnect(cid)
            canvas.mpl_disconnect(cidk)
    return peak_times,peak_vals

def onarrow(e):
    global scroll_slider
    if e.key == "right":
        pos = ax.get_xlim()
        scroll_slider.set_val(pos[0]+.500)
    elif e.key == "left":
        pos = ax.get_xlim()
        scroll_slider.set_val(pos[0]-.500)
    else:
        return

cid = fig.canvas.mpl_connect('button_press_event', onclick)
cidk = fig.canvas.mpl_connect('key_press_event', onarrow)







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

