import numpy as np
import tkinter as tk
from tkinter.filedialog import askopenfilename
import pandas as pd
from scipy.signal import find_peaks as findpeaks
from scipy.signal import firwin, lfilter
import bioread
import os
import matplotlib.pyplot as plt

def find_preceed_peak(signal,times,specific_time,preceed_window):
    idx = np.searchsorted(times, specific_time, side="left")
    local_max = np.argmax(signal[idx-np.round(preceed_window).astype(int):idx+1])
    global_max = idx-np.round(preceed_window).astype(int)+local_max
    return global_max,np.max(signal[idx-np.round(preceed_window).astype(int):idx+1])

def mm_lopass(data,times,srate,cutoff):
    tb = np.minimum(np.maximum(0.25 * cutoff, 2.),srate / 2. - cutoff)
    N = 3.3 * 1 / tb
    N = max(int(np.ceil(N *1*srate)), 1)
    N += (N - 1) % 2
    cutoff/(srate/2)
    b = firwin(N, cutoff/(srate/2),window='hamming')
    x = lfilter(b, 1.0, data)
    d = 0.5 * (N-1) / srate
    x=x[N-1:]
    t=times[N-1:]-d   
    return t,x

def loadem():
    
    def load_acq():
        root = tk.Tk()
        root.withdraw()
        root.update()
        file_path = askopenfilename()
        root.destroy()
        acq_dataset=bioread.read_file(file_path)
        return acq_dataset,file_path
    
    def select_chan(data,instr,prepop):
        window = tk.Tk()
        window.title(instr)
        window.geometry('300x600')
        
        var=[]
        c=[]
        for ind in range(len(data.channel_headers)):
            if ind==prepop[0]:
                var.append(tk.IntVar(value=1))
            else:
                var.append(tk.IntVar())
            c.append(tk.Checkbutton(window, text=data.channel_headers[ind].name,variable=var[ind], onvalue=1, offvalue=0))
            c[ind].pack()
        
        var_highp = tk.IntVar(value=prepop[1])
        c_highp = tk.Checkbutton(window, text="perform highpass filter",variable=var_highp, onvalue=1, offvalue=0)
        c_highp.pack(pady=10)
        
        var_deriv = tk.IntVar(value=prepop[2])
        c_deriv = tk.Checkbutton(window, text="compute a derivative",variable=var_deriv, onvalue=1, offvalue=0)
        c_deriv.pack(pady=20)
        
        exit_button = tk.Button(window, text="Continue...", command=window.destroy)
        exit_button.pack(pady=20)
        window.mainloop()
        chan=np.argwhere([v.get()==1 for v in var])
        highp=var_highp.get()==1
        deriv=var_deriv.get()==1
        return chan[0][0],deriv,highp
    
    def get_cont_ts(acq_dataset,ind,deriv,highp):
        hz=acq_dataset.channels[ind].samples_per_second
        t=acq_dataset.channels[ind].time_index
        s=acq_dataset.channels[ind].data
        
        if deriv==True:
            s=np.diff(s)
            t=t[1:]
            
        if highp==True:
            t,s=mm_lopass(s,t,hz,22.5)
            print('applying lowpass filter, 22.5 Hz cutoff')
        
        out_dict={}
        out_dict['s']=s
        out_dict['t']=t
        out_dict['hz']=hz
        
        return out_dict
    
    acq_dataset,file_path=load_acq()
    ind,deriv,highp=select_chan(acq_dataset,'Select Acceleration Channel to estimate time intervals',[3,1,0])
    print('selected contractility channel is '+acq_dataset.channel_headers[ind].name)
    print('perform highpass: '+str(highp))
    print('compute derivative: '+str(deriv))
    cont_dict=get_cont_ts(acq_dataset,ind,deriv,highp)
    
    return cont_dict,file_path
    
def compute_peaks(cont_dict):
    t=cont_dict['t']
    s=cont_dict['s']
    hz=cont_dict['hz']
    thresh_window = tk.Tk()
    thresh_window.withdraw()
    thresh_window.update()
    thresh = float(tk.simpledialog.askstring("Input", "Select max peak threshold",parent=thresh_window))
    thresh_window.destroy()
    p_ind,p_dict = findpeaks(s,thresh,distance=.5*hz)
    peak_times=t[p_ind]
    peak_vals=p_dict['peak_heights']
    return peak_times,peak_vals

def thresh_ind(hz):
    ind=np.arange(0,20*hz).astype(int)
    return ind

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

def shift_peaks(y,peak_times,peak_plot,cont_s,cont_t):
    old_peak_times=peak_times.copy()
    cont_vals=np.zeros(len(peak_times))
    peak_times=np.zeros(len(old_peak_times))
    for peak_ind,peak in enumerate(old_peak_times):
        new_peak_idx,cont_vals[peak_ind]=find_preceed_peak(cont_s,cont_t,peak,y)
        peak_times[peak_ind]=cont_t[new_peak_idx]
    peak_plot.set_xdata(peak_times)
    peak_plot.set_ydata(cont_vals)
    return peak_times,cont_vals,peak_plot

#def find_cont_from_acc(acc_peak_times,acc_timeseries,cont_time
#
# def confirm_peak_times(peak_times,peak_vals,cont_dict):
#   
#     def compute_hr(peak_times,peak_vals,cont_dict,window):
#         print('happy')
#        
#     confirm_window = tk.Tk()
#     confirm_window.title("Are you happy with time intervals?")
#     confirm_window.withdraw()
#     confirm_window.update()
#     yes_button = tk.Button(window, text="Yes", command=compute_hr(peak_times,peak_vals,cont_dict,window))
#     yes_button.pack()
#     no_button = tk.Button(window, text="No", command=window.destroy)
#     no_button.pack()
#     confirm_window.mainloop()