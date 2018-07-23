import numpy as np
import pandas as pd
import sys
import os
import Qdata
import pickle
import subprocess
import find_O2, pickle_qout

def catalyst_name(str): # returns name of bare catalyst given catalystO2 .out file name
    end = str.find("_") # we expect file name to contain "_" before job settings
    beg = end
    index = -1
    while(index == -1):
        if beg == -1:
            return None
        index = str.find("O2", beg, end)
        beg -= 1
    return str[0:index]

cat_dir, catO2_dir = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/')
# directories containing .out files

catO2_fn = catO2_dir.split('.')[0].split('/')[-1]
cat_fn = cat_dir.split('.')[0].split('/')[-1]

pickle_qout.pickle_data(catO2_dir, catO2_fn + ".p")
find_O2.create_df(catO2_fn + ".p")
catO2_df = pickle.load(open(catO2_fn + "_df.p", "rb"))
pickle_qout.pickle_data(cat_dir, cat_fn + ".p")
cat_qdata = pickle.load(open(cat_fn + ".p", "rb"))

cat_list = []
for entry in catO2_df["File_Name"]:
    cat_list.append(catalyst_name(entry))
cat_fn_list, cat_energy_list = [], []
for entry in cat_list:
    for qdata in cat_qdata:
        if qdata.filename.split('_')[0] == entry:
            cat_fn_list.append(qdata.filename)
            cat_energy_list.append(qdata.E)
            break
            # assumes cat_dir does not contain multiple files for the same catalyst
    else:
        cat_fn_list.append(None)
        cat_energy_list.append(None)

cat_df = pd.DataFrame.from_items([('Catalyst_File_Name', cat_fn_list), ('Catalyst_Energy', cat_energy_list)])
matched_df = pd.concat([cat_df, catO2_df], axis = 1)
matched_df.rename(index = str, columns = {"File_Name" : "CatalystO2_File_Name", "Energy" : "CatalystO2_Energy"}, inplace = True)

matched_df.to_csv(cat_fn + '_matched.csv')
matched_df.to_pickle(cat_fn + '_matched_df.p')
