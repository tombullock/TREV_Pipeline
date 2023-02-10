TUTORIAL

Preprocessing the example dataset was tested with Jupyter notebook version 3.8.12, installed via Anaconda version 4.10.3
No further installation of packages should be necessary*

Example dataset can be downloaded here: https://ucsb.box.com/s/qewjv8nlq8kaoxicvp75h96ssv3dmrs1

To run the preprocessing script:

- Clone the TREV_Pipeline github repo: https://github.com/caitgregory/TREV_Pipeline
- Open miniMEAP.ipynb in Jupyter notebook
- Run all 5 cells:

Cell 1 loads the data via a GUI:
- user locates the AcqKnolwedge file (example dataset: IV_301_1.acq)
- user specifies which is acceleration channel (at the top) and which is respiration channel (at the bottom)
- user specfiies if respiration channel needs lowpass (if using z0; not needed if you have the resp belt).
- click "continue", which closes the GUI
- give this a minute to load/preprocess everything, "data loaded!" will appear in the progress bar
- notebook currently preset for example dataset, i.e., channels for acceleration and respiration will already be selected, along with option to lowpass respiration. 

Cell 2: 
- plots a section (20 seconds worth) of acceleration peaks to pick a threshold for min amplitude.
- enter in the accompanying GUI and click OK (will close the plot)
- (in the example dataset, I use 0.40. There are spurious smaller peaks, but they occur <.5 s after each larger peak (possible P.A.C.?)
- Err on a lower threshold

Cell 3:
- plots an interactive timeseries of acceleration. Note this is just to identify beat times. Acceleration peaks are more robust to noise.
- You can remove, add and estimate (moving ensemble average) peaks with:
- left click - add peak
- right click - remove peak
- hold m and left click - use moving average to estimate peak and interval (remove any unwanted peak first)
- Close the window manually when you're done.

Cell 4:
- plots an interactive timeseries of contractility (derivative of acceleration).
- You can again remove, add and estimate (moving ensemble average) peaks with:
- left click - add peak
- right click - remove peak
- hold m and left click - use moving average to estimate peak and interval
- Close the window manually when you're done.

Cell 5: 
- removes respiration and heart rate from data and saves a csv of a timestamp and contractility at each heartbeat
- example dataset I get 1093 beats

*Cell 1 installs bioread via pip. pip was preloaded in my version of Jupyter notebook, but you might need to install separately....