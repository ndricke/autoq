import numpy as np
import pandas as pd
import sys
import os
import Qdata
import pickle
import subprocess
import find_O2, pickle_qout
from collections import OrderedDict

def catalyst_name(str): # returns name of bare catalyst given catalystO2 .out file name
    str_bn = os.path.basename(str)
    return str_bn[0:str_bn[0:str_bn.find("_")].rfind("O2")]

def dataPull(cat_list, cat_qdata):
    cat_fn_list, cat_energy_list, cat_AS_chelpg_list = [], [], []
    for ind, entry in enumerate(cat_list):
        for qdata in cat_qdata:
            if qdata.filename.split('_')[0] == entry:
                cat_fn_list.append(qdata.filename)
                try:
                    cat_energy_list.append(qdata.E)
                except:
                    cat_energy_list.append(None)
                try:
                    chelpgs = [float(item) for item in qdata.chelpg]
                    active_site_index = catO2_df.iloc[ind]['Active_Site']
                    cat_AS_chelpg_list.append(chelpgs[active_site_index - 1])
                except:
                    cat_AS_chelpg_list.append(None)
                break
                # assumes cat_dir does not contain multiple files for the same catalyst
        else:
            cat_fn_list.append(None)
            cat_energy_list.append(None)
            cat_AS_chelpg_list.append(None)

    return cat_fn_list, cat_energy_list, cat_AS_chelpg_list


## The process of pickling data should only be done if we can't find the appropriate pickle files
# directories containing .out files
cat_dir, catO2_dir, cat_c1_dir = sys.argv[1].rstrip('/'), sys.argv[2].rstrip('/'), sys.argv[3].rstrip('/')

catO2_fn = catO2_dir.split('.')[0].split('/')[-1]
cat_fn = cat_dir.split('.')[0].split('/')[-1]
cat_c1_fn = cat_c1_dir.split('.')[0].split('/')[-1]



##Generate and pickle a list of Qdata objects from parsing all files in catO2_dir and cat_dir
pickle_qout.pickle_data(catO2_dir, catO2_fn + ".p")
pickle_qout.pickle_data(cat_dir, cat_fn + ".p")
pickle_qout.pickle_data(cat_c1_dir, cat_c1_fn + ".p")

find_O2.create_df(catO2_fn + ".p") #use find_O2 to generate a dataframe from a list of Qdata objects


##Load data if previously generated
catO2_df = pickle.load(open(catO2_fn + "_df.p", "rb"))
cat_qdata = pickle.load(open(cat_fn + ".p", "rb"))
cat_c1_qdata = pickle.load(open(cat_c1_fn + ".p", "rb"))

cat_list = []
for entry in catO2_df["CatalystO2_File_Name"]:
    cat_list.append(catalyst_name(entry))

# collect a0 bare catalyst data
cat_fn_list, cat_energy_list, cat_AS_chelpg_list = [], [], []
for ind, entry in enumerate(cat_list):
    for qdata in cat_qdata:
        if qdata.filename.split('_')[0] == entry:
            cat_fn_list.append(qdata.filename)
            try:
                cat_energy_list.append(qdata.E)
            except:
                cat_energy_list.append(None)
            try:
                chelpgs = [float(item) for item in qdata.chelpg]
                active_site_index = catO2_df.iloc[ind]['Active_Site']
                cat_AS_chelpg_list.append(chelpgs[active_site_index - 1])
            except:
                cat_AS_chelpg_list.append(None)
            break
            # assumes cat_dir does not contain multiple files for the same catalyst
    else:
        cat_fn_list.append(None)
        cat_energy_list.append(None)
        cat_AS_chelpg_list.append(None)

# collect c1 bare catalyst data
cat_c1_fn_list, cat_c1_energy_list, cat_c1_AS_chelpg_list = [], [], []
for ind, entry in enumerate(cat_list):
    for qdata in cat_c1_qdata:
        if qdata.filename.split('_')[0] == entry:
            cat_c1_fn_list.append(qdata.filename)
            try:
                cat_c1_energy_list.append(qdata.E)
            except:
                cat_c1_energy_list.append(None)
            try:
                chelpgs = [float(item) for item in qdata.chelpg]
                active_site_index = catO2_df.iloc[ind]['Active_Site']
                cat_c1_AS_chelpg_list.append(chelpgs[active_site_index - 1])
            except:
                cat_c1_AS_chelpg_list.append(None)
            break
    else:
        cat_c1_fn_list.append(None)
        cat_c1_energy_list.append(None)
        cat_c1_AS_chelpg_list.append(None)

cat_dict = {'Catalyst_File_Name': cat_fn_list, 'Catalyst_Energy': cat_energy_list, 'Catalyst_Active_Site_CHELPG': cat_AS_chelpg_list}
cat_df = pd.DataFrame.from_dict(OrderedDict(cat_dict))
cat_c1_dict = {'Catalyst_c1_File_Name': cat_c1_fn_list, 'Catalyst_c1_Energy': cat_c1_energy_list, 'Catalyst_c1_Active_Site_CHELPG': cat_c1_AS_chelpg_list}
cat_c1_df = pd.DataFrame.from_dict(OrderedDict(cat_c1_dict))
matched_df = pd.concat([cat_df, catO2_df, cat_c1_df], axis = 1)

matched_df.to_csv(cat_fn + '_matched.csv')
matched_df.to_pickle(cat_fn + '_matched_df.p')
