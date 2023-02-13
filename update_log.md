update 12/7/22 7:40 - nmd
- now loads in from .acq files
- now uses a gui to allow user select contractility channel
- now uses a gui to allow user select respiration channel
- cell 2 now plots 10 s of data to let user specify amp threshold in findpeaks
- minor tweaks to x and y zoom on interactive plots

update 12/8/22 10:30 - nmd
- minor tweaks here and there - considering using dZ for temp, then those temp indices for cont amps

update 2/9/22 6:48 - nmd
- beta version should be good to go
- using acceleration channel to identify heartbeats (equivalent R intervals)
- Cell 1 loads the data via a GUI
- user specifies which is acceleration channel and which is respiration channel
- user specfiies if respiration channel needs lowpass (if using z0)
- give this a minute to load/preprocess everything, "data loaded!" will appear in progress bar
- first plot (Cell 2) plots a section (20 seconds worth) of acceleration peaks to pick a threshold for min amplitude. Enter in the accompanying GUI.
- second plot (Cell 3) interactive plot of entire series of acceleration peaks to adjust errors in peak identification (see instructions in cell)
- third plot (Cell 4) replots contractility and allows a final adjustment
- Cell 5 models out respiration (volume and phase) and saves a csv. A column of peak times and contractility.

- adding tutorial.md
- deleting various old files

update 2/10/22 -nmd
- change to miniMEAPsys.py, can now work with a respiration belt, if specified
- change to miniMEAPsys.py, user now specifies an output filename, prepopulated with default

update 2/13/22 - nmd
- adding new function for cell 4, adjust_peak_amp, that only changes peak amplitude
- adding new function for cell 4, meap_dat_cell4, that only changes peak amplitude
- Cell 3: plot now updates number of peaks as they are added and removed
- loadem now does a plt.closeall and a 
- bug fix - find_preceed_peak - now searches either in the preceeding 250 datapoints for cont peak, or if peak ind <250, the entire series up to that point
- bug fix - meap_dat - now boolean test to make sure there are two peaks either side of meap_dat loc
- bug fix - duplicate shift_peaks removed

- desired future fixes: 
some kind of tkinter "kill all / close all" that can be called in Cell 1, might help with loading issues
some way of making functions of the text in Cells 3 and 4, which make the notebook very long
Cell 1 very buggy, need to restart kernel a lot... 

