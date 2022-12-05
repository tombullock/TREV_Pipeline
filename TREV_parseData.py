
# TREV_parseData.py
# Caitlin Gregory, UCSB Attention Lab
# Nov 21st, 2022

# Purpose:
#   load in Acqknowledge files in .txt format,
#   parse data with respect to experimental triggers,
#   and save epoched data into separate files

# ---

import pandas as pd
import numpy as np
import os

sj_num = input('enter the sj number you want to process:   ') # a string

# set up directories
parent_dir = '/Users/caitlingregory/Desktop/Projects/RADARS/python/'
data_dir = parent_dir + 'data_raw/'
dest_dir = parent_dir + 'data_parsed/'
this_sj_dir = dest_dir + 'sj' + sj_num
filename = data_dir + 'sj' + sj_num + '_RADARS.txt'

# load in data (raw data in .txt, output into pd dataframe)
data = pd.read_csv(filename, delimiter=',', header=None)

# name the columns (this is hard coded based on RADARS template)
data.columns = ['z', 'dz', 'resp_belt', 'bp',
                'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                'z_filtered', 'dZ_filtered', 'exp_trig', 'acc', 'resp_TREV']

# drop the columns named 'NA'
data.drop(columns=[col_name for col_name in data.columns if 'NA' in col_name], inplace=True)

# check length of data, should be around 90-120 minutes (radars)
t = np.arange(0, len(data)) / 1000

# epoch the data now
#   goal: parse out breaks and between block times
#   then separate data depending on block type

data.reset_index(inplace=True, drop=True)

# quick reality check. only continue if session length looks correct
print('session appears to be ', round((t[-1] / 60), 2), ' minutes long...')
good = input('continue?  ')
if good == 'no':
    print('please reload the data and start again')
elif good == 'yes':

    # epoch the data

    all_triggers_idx = data.index[data['exp_trig'] != 0].tolist()  # get indexes of each trigger
    data_orig = data.copy();  # make copy of data just incase 
    os.makedirs(this_sj_dir, exist_ok=True)

    # list triggers
    pre_post_trigs = [58, 61, 62, 63, 64, 65, 81, 82, 83, 84, 85]
    cpt_trigs = [11, 12, 13, 21, 22, 23, 31, 32, 33, 41, 42, 43, 51, 52, 53]

    # get idx of trig64 (there was a trigger error in radars with too many 64s)
    trig64_idx = data[data['exp_trig'] == 64].index

    for this_idx in all_triggers_idx:

        this_code = data['exp_trig'][this_idx]
        this_file = '/sj' + sj_num + '_bl' + str(this_code) + '_physio.csv'

        
        if this_code in pre_post_trigs and this_idx != trig64_idx[0]: # second conditional radars specific

            if os.path.exists(this_sj_dir + this_file) == False:
                this_dur = 210000  # 210s
                data_new = data.loc[range(this_idx, this_idx + this_dur)]

                # save block data into sj dir
                print('parsing block ' + str(this_code))
                data_new.to_csv(this_sj_dir + this_file)

        elif this_code in cpt_trigs:

            if os.path.exists(this_sj_dir + this_file) == False:
                this_dur = 90000  # 90s
                data_new = data.loc[range(this_idx, this_idx + this_dur)]

                # save block data into sj dir
                print('parsing block ' + str(this_code))
                data_new.to_csv(this_sj_dir + this_file)

print('done')

