import numpy as np
import pandas as pd
import sys
import os
import pickle
import subprocess
from collections import OrderedDict
import argparse

import Qdata
import find_O2, pickle_qout

def catalyst_name(filename): # returns name of bare catalyst given catalystO2 .out file name
    end = filename.find("_") # we expect file name to contain "_" before job settings
    beg = end
    index = -1
    while(index == -1):
        if beg == -1:
            return None
        index = filename.find("O2", beg, end)
        beg -= 1
    return filename[0:index]

def GenLoadCat(cat_dir, save=False):
    qdata_list = pickle_qout.collect_qdata(cat_dir)
    if save == True:
        cat_fn = cat_dir.split('.')[0].split('/')[-1]
        pickle.dump(qdata_list, open(cat_fn+'.p', 'wb'))
    return qdata_list

def GenLoadCatO2(catO2_dir, save=False):
    catO2_fn = catO2_dir.split('.')[0].split('/')[-1]
    pickle_qout.pickle_data(catO2_dir, catO2_fn + ".p")
    O2_df = find_O2.create_df(catO2_fn + ".p", save_df=save) #use find_O2 to generate a dataframe from a list of Qdata objects
    return O2_df

#Could we directly add the data from each qdata to the dataframe rather than adding to a list?    
## What needs to happen is, for an entry in the catO2_df, find the corresponding bare catalysts qdata
class MatchO2(object):

    def __init__(self, in_O2_df, cat_qdata, cation_qdata):
    
        energy = "Catalyst_Energy"
        cat_fn = "Catalyst_File_Name"
        AS_CHELPG = "Catalyst_Active_Site_CHELPG"
        c1_energy = "Catalyst_c1_Energy"
        c1_cat_fn = "Catalyst_c1_File_Name"
        c1_AS_CHELPG = "Catalyst_c1_Active_Site_CHELPG"
    
        self.catO2_df = in_O2_df.copy()
        self.catO2_df = self.catO2_df.assign(Catalyst_File_Name=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_Energy=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_Active_Site_CHELPG=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_c1_File_Name=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_c1_Energy=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_c1_Active_Site_CHELPG=None)

        self.FillMatchedO2(cat_qdata, energy, cat_fn, AS_CHELPG)
        self.FillMatchedO2(cation_qdata, c1_energy, c1_cat_fn, c1_AS_CHELPG)


    def FillMatchedO2(self, qdatas, energy, cat_fn, AS_CHELPG):
        for i, entry in enumerate(self.catO2_df["CatalystO2_File_Name"]):
            catO2_entry_name = catalyst_name(entry)
            for qdata in qdatas:
                #if qdata.filename.split('_')[0] == catO2_entry_name:
                if qdata.filename.split('_')[1] == entry.split('_')[1]:
                    self.catO2_df.at[i, energy] = qdata.E
                    self.catO2_df.at[i, cat_fn] = qdata.filename
                    try:
                        chelpgs = [float(item) for item in qdata.chelpg]
                        active_site_index = self.catO2_df.iloc[i]['Active_Site']
                        self.catO2_df.at[i,AS_CHELPG] = chelpgs[active_site_index - 1]
                    except:
                        print("Could not get CHELPG charges for " + qdata.filename)
                    break # assumes cat_dir does not contain multiple files for the same catalyst


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-bare", help="Bare catalyst in active charge/spin state", type=str)
    parser.add_argument("-O2", help="Catalyst with O2 bound", type=str)
    parser.add_argument("-cation", help="Ionized bare catalyst", type=str)
    parser.add_argument("-save", help="Save intermediate qdata_list pickles", type=bool, default=False)
    args = parser.parse_args()

    if (os.path.isfile(args.bare+'.p') and os.path.isfile(args.cation+'.p') and os.path.isfile(args.O2+'_df.p')):
        cat_qdata_list = pickle.load(open(args.bare + ".p", "rb"))
        cation_qdata_list = pickle.load(open(args.cation + ".p", "rb"))
        O2_df = pickle.load(open(args.O2 + "_df.p", "rb"))
    else:
        cat_qdata_list = GenLoadCat(args.bare, save=args.save)
        cation_qdata_list = GenLoadCat(args.cation, save=args.save)
        O2_df = GenLoadCatO2(args.O2, save=args.save)

    match = MatchO2(O2_df, cat_qdata_list, cation_qdata_list)
    O2_matched_df = match.catO2_df
    print(O2_matched_df)

    O2_matched_df.to_csv('catO2_matched.csv')



        













    

