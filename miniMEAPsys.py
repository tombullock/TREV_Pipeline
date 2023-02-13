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

def outsave_box(file_in):
    out_window = tk.Tk()
    out_window.geometry('400x300+500+500')
    out_window.withdraw()
    out_window.update()
    out_f = tk.simpledialog.askstring("Save File", "Select output filename",parent=out_window,initialvalue=file_in)
    out_window.destroy()
    return out_f

def regress_out(cont_dict):
    nb = len(cont_dict['cont_peak_times'])
    RR = cont_dict['RR']
    resp_amount = 1000*np.ones(nb)
    resp_cycle = 1000*np.ones(nb)
    for beat_no,beat_time in enumerate(cont_dict['cont_peak_times']):    
        idx = np.searchsorted(cont_dict['resp_t'], beat_time, side="left")
        if idx!=len(cont_dict['resp_t']):
            if np.abs(beat_time-cont_dict['resp_t'][idx])<.01:
                resp_amount[beat_no]=cont_dict['resp_amount'][idx]
                resp_cycle[beat_no]=cont_dict['resp_cycle'][idx]

    resp_amount[resp_amount==1000]=np.median(resp_amount)
    resp_cycle[resp_cycle==1000]=np.median(resp_cycle)
    RR[np.isnan(RR)]=np.nanmedian(RR)
    X=np.array([np.ones(nb),zscore(RR),zscore(resp_amount),zscore(resp_cycle)]).T
    y = cont_dict['raw_contractility_peaks']
    B = np.matmul(np.linalg.inv(np.matmul(X.T,X)),np.matmul(X.T,y))
    resid = y - np.matmul(X,B)
    cont_dict['resid_contractility'] = B[0]+resid
    return cont_dict

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y

def find_preceed_peak(signal,times,specific_time,preceed_window):
    idx = np.searchsorted(times, specific_time, side="left")
    if idx>preceed_window:
        local_max = np.argmax(signal[idx-np.round(preceed_window).astype(int):idx+1])
        global_max = (idx-np.round(preceed_window).astype(int)+local_max).astype(int)
        amp = np.max(signal[idx-np.round(preceed_window).astype(int):idx+1])
    else:
        local_max = np.argmax(signal[0:idx+1])
        global_max = local_max
        amp = np.max(signal[0:idx+1])
    return global_max,amp

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

def norm_range(data,new_min,new_max):
    new_range = new_max-new_min
    Z = (data-np.min(data))/(np.max(data)-np.min(data))
    norm_data = new_min + Z*new_range
    return norm_data

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
        window.geometry('300x800')
        
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
        c_highp = tk.Checkbutton(window, text="estimate resp with low-pass filter",variable=var_highp, onvalue=1, offvalue=0)
        c_highp.pack(pady=20)
        
        rvar=[]
        rc=[]
        for ind in range(len(data.channel_headers)):
            if ind==prepop[2]:
                rvar.append(tk.IntVar(value=1))
            else:
                rvar.append(tk.IntVar())
            rc.append(tk.Checkbutton(window, text=data.channel_headers[ind].name,variable=rvar[ind], onvalue=1, offvalue=0))
            rc[ind].pack()
               
        exit_button = tk.Button(window, text="Continue...", command=window.destroy)
        exit_button.pack(pady=20)
        window.mainloop()
        acc_chan=np.argwhere([v.get()==1 for v in var])
        resp_chan=np.argwhere([v.get()==1 for v in rvar])
        highp=var_highp.get()==1
        return acc_chan[0][0],resp_chan[0][0],highp
    
    def get_cont_ts(acq_dataset,acc_chan,resp_chan,highp):
        ##Contractility
        hz=acq_dataset.channels[acc_chan].samples_per_second
        t=acq_dataset.channels[acc_chan].time_index
        s=acq_dataset.channels[acc_chan].data
        print('removing slow movements from acceleration data with 7-order polynomial')
        s_p = np.polyfit(t,s,10)
        s = s-np.polyval(s_p,t)
        print('applying lowpass filter to acceleration channel, 22.5 Hz cutoff')
        t,s=mm_lopass(s,t,hz,22.5)       
        
        resp_hz=acq_dataset.channels[resp_chan].samples_per_second
        resp_t=acq_dataset.channels[resp_chan].time_index
        resp_s=acq_dataset.channels[resp_chan].data
        print('removing slow movements from respiration data with 7-order polynomial')
        resp_s_p = np.polyfit(resp_t,resp_s,10)
        resp_s = resp_s-np.polyval(resp_s_p,resp_t)
        if highp==1:
            print('applying lowpass filter to respiration channel, 0.35 Hz cutoff')
            resp_t,resp_s=mm_lopass(resp_s,resp_t,resp_hz,0.35)
        else:
            print('NOT applying lowpass to respiration channel')
                        
        out_dict={}
        out_dict['s']=s
        out_dict['t']=t
        out_dict['hz']=hz
        
        out_dict['raw_resp_s']=resp_s
        out_dict['raw_resp_t']=resp_t
        out_dict['raw_resp_hz']=resp_hz
        
        return out_dict
    
    def resp_cycle_amount(cont_dict, min_dist):
        peak_inds,_ = findpeaks(cont_dict['raw_resp_s'],distance=min_dist)
        k=1
        peak_no=-1
        resp_t=[]
        resp_cycle=[]
        resp_amount=[]
        while k>0:
            peak_no=peak_no+1
            if peak_no+2<len(peak_inds):
                next3=cont_dict['raw_resp_s'][peak_inds[np.arange(peak_no,peak_no+3)]]

                if (next3[0]<next3[1]) & (next3[2]<next3[1]):
                    wave_t = cont_dict['raw_resp_t'][peak_inds[peak_no]:peak_inds[peak_no+2]]
                    resp_t.append(wave_t)
                
                    up_wave = cont_dict['raw_resp_t'][peak_inds[peak_no]:peak_inds[peak_no+1]]
                    up_phase = -np.sin(norm_range(up_wave,0,np.pi))
                
                    dn_wave = cont_dict['raw_resp_t'][peak_inds[peak_no+1]:peak_inds[peak_no+2]]
                    dn_phase = -np.sin(norm_range(dn_wave,np.pi,2*np.pi))
                
                    phase = np.concatenate((up_phase,dn_phase))
                    resp_cycle.append(phase)
                
                    wave_vals = cont_dict['raw_resp_s'][peak_inds[peak_no]:peak_inds[peak_no+2]]
                    resp_amount.append(wave_vals-np.mean(wave_vals))
                
                    peak_no=peak_no+1
            else:
                k=0
    
        cont_dict['resp_t']=np.hstack(resp_t)
        cont_dict['resp_amount']=np.hstack(resp_amount)
        cont_dict['resp_cycle']=np.hstack(resp_cycle)
        cont_dict['resp_peak_inds']=peak_inds
    
        return cont_dict
    plt.close('all')  
    acq_dataset,file_path=load_acq()
    acc_chan,resp_chan, highp=select_chan(acq_dataset,'Select Acceleration and Respiration Channels',[3,1,1])
    print('selected contractility channel is '+acq_dataset.channel_headers[acc_chan].name)
    print('selected respiration channel is '+acq_dataset.channel_headers[resp_chan].name)
    
    print('extracting raw signals from AcqKnowledge files...')
    cont_dict=get_cont_ts(acq_dataset,acc_chan,resp_chan,highp)
    
    print('estimating respiration amount and cycle...')
    cont_dict=resp_cycle_amount(cont_dict,cont_dict['raw_resp_hz']*.75)
    
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
    return peak_times,peak_vals,p_ind

def thresh_ind(hz):
    ind=np.arange(0,20*hz).astype(int)
    return ind

def meap_dat(ix,iy,peak_times,peak_vals,peak_plot):
    ind_after=np.argmax(peak_times>ix)
    if (ind_after-2>-1) | (ind_after+2<=len(peak_times)):
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

def meap_dat_cell4(ix,iy,peak_times,peak_vals,peak_plot):
    idx = np.argmin(np.abs(ix-peak_times))
    if (idx-2>-1) | (idx+2<=len(peak_times)):
        mu_inds = np.arange(idx-2,idx+2)
        new_val = np.mean(peak_vals[mu_inds])
        peak_vals[idx] = new_val
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

def adjust_peak_amp(ix,iy,peak_times,peak_vals,peak_plot):
    idx = np.argmin(np.abs(ix-peak_times))
    peak_vals[idx] = iy
    peak_plot.set_xdata(peak_times)
    peak_plot.set_ydata(peak_vals)
    return peak_times,peak_vals,peak_plot

def shift_peaks(yzoom,peak_times,cont_t,cont_s,cont_inds,peak_plot):
    if yzoom>=1:
        new_cont_inds=np.zeros(len(peak_times)).astype(int)
        for peak_ind,peak in enumerate(peak_times):
            new_cont_inds[peak_ind],_=find_preceed_peak(cont_s,cont_t,peak,yzoom)
    else:
        new_cont_inds=peak_inds.copy()
    peak_plot.set_xdata(cont_t[new_cont_inds])
    peak_plot.set_ydata(cont_s[new_cont_inds])
    ind_offset = cont_inds-new_cont_inds
    return ind_offset

def shift_peaks_no_plot(yzoom,peak_times,cont_t,cont_s,cont_inds):
    new_cont_inds=np.zeros(len(peak_times)).astype(int)
    for peak_ind,peak in enumerate(peak_times):
        new_cont_inds[peak_ind],_=find_preceed_peak(cont_s,cont_t,peak,yzoom)

    return new_cont_inds

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