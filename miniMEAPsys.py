import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import pandas as pd
from scipy.signal import find_peaks as findpeaks

def loadem():
    Tk().withdraw() 
    in_filename = askopenfilename() 
    df = pd.read_csv(in_filename,delimiter='\t',header=None)
    df.columns=['tom','z','dz','resp_belt','bp','tom','tom','tom','tom',\
    'tom','tom','tom','tom','z_filtered','dZ_filtered','exp_trig','acc','resp_TREV','tom']
    df.drop(columns=[s for s in df.columns if 'tom' in s],inplace=True)
    file_state=0
    if 'df' in locals():
        file_state=in_filename+' successfully loaded'
    return df,file_state

def meap_dat(ix,iy,peak_times,peak_vals,peak_plot):
    ind_after=np.argmax(peak_times>ix)
    mu_inds = np.arange(ind_after-2,ind_after+2)
    new_time0 = peak_times[mu_inds[1]]+np.diff(peak_times[mu_inds[0:2]])
    new_time1 = peak_times[ind_after]-np.diff(peak_times[mu_inds[2:]])
    new_time = np.mean([new_time0,new_time1])
    new_val = np.mean(peak_vals[mu_inds])
    peak_times = np.insert(peak_times, ind_after-1, new_time)
    peak_vals = np.insert(peak_vals, ind_after-1, new_val)
    peak_plot.set_xdata(peak_times)
    peak_plot.set_ydata(peak_vals)
    return peak_times,peak_vals,peak_plot

def remove_point(ix,peak_times,peak_vals,peak_plot):
    ind=np.argmin(np.abs(ix-peak_times))
    peak_times = np.delete(peak_times,ind)
    peak_vals = np.delete(peak_vals,ind)
    peak_plot.set_xdata(peak_times)
    peak_plot.set_ydata(peak_vals)
    return peak_times,peak_vals,peak_plot
    
def add_point(ix,iy,peak_times,peak_vals,peak_plot):
    ind=np.argmax(peak_times>ix)
    peak_times = np.insert(peak_times, ind, ix)
    peak_vals = np.insert(peak_vals, ind, iy)
    peak_plot.set_xdata(peak_times)
    peak_plot.set_ydata(peak_vals)
    return peak_times,peak_vals,peak_plot

def get_ts(df):
    #if 'df' in globals():
    s = np.diff(df.acc)
    t = np.array(df.index[1:])/1000
    #else:
        #t = np.linspace(0.0, 100.0, 10000)
        #s = np.random.normal(size=10000)
    p_ind,p_dict = findpeaks(s,0.05,distance=.400)
    peak_times=t[p_ind]
    peak_vals=p_dict['peak_heights']
    return t,s,peak_times,peak_vals