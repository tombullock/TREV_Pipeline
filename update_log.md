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