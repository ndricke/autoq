import numpy as np
import pandas as pd
import sys
import os
import dill
import subprocess
from collections import OrderedDict
import argparse
import re

import Qdata
import find_O2, dill_qout

# check the file name format if this doesn't work
def catalyst_name(filename): # returns name of bare catalyst given catalystO2 .out file name
    return filename[0:filename[0:filename.find("_")].rfind("O2")]

# to match files of same catalyst with different bound species
# run name_conv.py first if necessary
def parse_filename(filename):
    if True: # janky workaround for SciFinder file names; make False otherwise
        if "-" in filename:
            return [filename[2:filename.find('x')], filename[filename.find('x') + 1], filename[filename.find('-') + 1:filename.find('_')]]
        else:
            return [filename[2:filename.find('x')], filename[filename.find('x') + 1], None]
    fn_fragments = filename.split('-')
    cat_type = fn_fragments[0]
    func_num = re.findall(r'\d+', fn_fragments[1])[0]
    if len(fn_fragments) == 3:
        active_site = re.findall(r'\d+', fn_fragments[2])[0]
    else:
        active_site = None
    return [cat_type, func_num, active_site] #strings

# at this point it probably doesn't matter but it might be more robust to make
    # parse_filename() find bound species type and then match active site only
    # if fn1 or fn2 is bare
def same_catalyst(fn1, fn2, match_active_site):
    fn1_parsed, fn2_parsed = parse_filename(fn1), parse_filename(fn2)
    if match_active_site:
        if fn1_parsed == fn2_parsed:
            return True
    else:
        if fn1_parsed[0] == fn2_parsed[0] and fn1_parsed[1] == fn2_parsed[1]:
            return True
    return False

def GenLoadCat(cat_dir, save=False):
    if cat_dir == None:
        return None
    qdata_list = dill_qout.collect_qdata(cat_dir)
    if save == True:
        cat_fn = cat_dir.rstrip('/').split('/')[-1].split('.')[0] # sometimes the . is a -? not standardized
        dill.dump(qdata_list, open(cat_fn+'.p', 'wb'))
    return qdata_list

def GenLoadCatO2(catO2_dir, save=False):
    catO2_fn = catO2_dir.rstrip('/').split('/')[-1].split('.')[0]
    dill_qout.pickle_data(catO2_dir, catO2_fn + ".p")
    O2_df = find_O2.create_df(catO2_fn + ".p", save_df=save) #use find_O2 to generate a dataframe from a list of Qdata objects
    return O2_df

#Could we directly add the data from each qdata to the dataframe rather than adding to a list?
## What needs to happen is, for an entry in the catO2_df, find the corresponding bare catalysts qdata
class MatchO2(object):

    def __init__(self, in_O2_df, cat_qdata, cation_qdata, OOH_qdata, OH_qdata, O_qdata, CO_qdata) : #, CN_qdata):

        energy, cat_fn, AS_CHELPG = "Catalyst_Energy", "Catalyst_File_Name", "Catalyst_Active_Site_CHELPG"
        c1_energy, c1_cat_fn, c1_AS_CHELPG = "Catalyst_c1_Energy", "Catalyst_c1_File_Name", "Catalyst_c1_Active_Site_CHELPG"
        OOH_energy, OOH_cat_fn, OOH_AS_CHELPG = "CatalystOOH_Energy", "CatalystOOH_File_Name", "CatalystOOH_Active_Site_CHELPG"
        O_energy, O_cat_fn, O_AS_CHELPG = "CatalystO_Energy", "CatalystO_File_Name", "CatalystO_Active_Site_CHELPG"
        OH_energy, OH_cat_fn, OH_AS_CHELPG = "CatalystOH_Energy", "CatalystOH_File_Name", "CatalystOH_Active_Site_CHELPG"
        CO_energy, CO_cat_fn, CO_AS_CHELPG = "CatalystCO_Energy", "CatalystCO_File_Name", "CatalystCO_Active_Site_CHELPG"
        #CN_energy, CN_cat_fn, CN_AS_CHELPG = "CatalystCN_Energy", "CatalystCN_File_Name", "CatalystCN_Active_Site_CHELPG"

        self.catO2_df = in_O2_df.copy()
        self.catO2_df = self.catO2_df.assign(Catalyst_File_Name=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_Energy=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_Active_Site_CHELPG=None)

        self.catO2_df = self.catO2_df.assign(Catalyst_c1_File_Name=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_c1_Energy=None)
        self.catO2_df = self.catO2_df.assign(Catalyst_c1_Active_Site_CHELPG=None)

        self.catO2_df = self.catO2_df.assign(CatalystOOH_File_Name=None)
        self.catO2_df = self.catO2_df.assign(CatalystOOH_Energy=None)
        self.catO2_df = self.catO2_df.assign(CatalystOOH_Active_Site_CHELPG=None)

        self.catO2_df = self.catO2_df.assign(CatalystO_File_Name=None)
        self.catO2_df = self.catO2_df.assign(CatalystO_Energy=None)
        self.catO2_df = self.catO2_df.assign(CatalystO_Active_Site_CHELPG=None)

        self.catO2_df = self.catO2_df.assign(CatalystOH_File_Name=None)
        self.catO2_df = self.catO2_df.assign(CatalystOH_Energy=None)
        self.catO2_df = self.catO2_df.assign(CatalystOH_Active_Site_CHELPG=None)

        self.catO2_df = self.catO2_df.assign(CatalystCO_File_Name=None)
        self.catO2_df = self.catO2_df.assign(CatalystCO_Energy=None)
        self.catO2_df = self.catO2_df.assign(CatalystCO_Active_Site_CHELPG=None)

        #self.catO2_df = self.catO2_df.assign(CatalystCN_File_Name=None)
        #self.catO2_df = self.catO2_df.assign(CatalystCN_Energy=None)
        #self.catO2_df = self.catO2_df.assign(CatalystCN_Active_Site_CHELPG=None)

        self.FillMatchedO2(cat_qdata, energy, cat_fn, AS_CHELPG, False)
        self.FillMatchedO2(cation_qdata, c1_energy, c1_cat_fn, c1_AS_CHELPG, False)
        self.FillMatchedO2(OOH_qdata, OOH_energy, OOH_cat_fn, OOH_AS_CHELPG, True)
        self.FillMatchedO2(O_qdata, O_energy, O_cat_fn, O_AS_CHELPG, True)
        self.FillMatchedO2(OH_qdata, OH_energy, OH_cat_fn, OH_AS_CHELPG, True)
        self.FillMatchedO2(CO_qdata, CO_energy, CO_cat_fn, CO_AS_CHELPG, True)
        #self.FillMatchedO2(CN_qdata, CN_energy, CN_cat_fn, CN_AS_CHELPG, True)

    def FillMatchedO2(self, qdatas, energy, cat_fn, AS_CHELPG, not_bare):
        if qdatas == None:
            return None # the DF columns will exist but be empty
        for i, entry in enumerate(self.catO2_df["CatalystO2_File_Name"]):
            for qdata in qdatas:
                if same_catalyst(entry, qdata.filename, not_bare):
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

    parser = argparse.ArgumentParser() # directories or .p files, I think
    parser.add_argument("-bare", help="Bare catalyst in active charge/spin state", type=str)
    parser.add_argument("-O2", help="Catalyst with O2 bound", type=str)
    parser.add_argument("-OOH", help="Catalyst with OOH bound", type=str)
    parser.add_argument("-OH", help="Catalyst with OH bound", type=str)
    parser.add_argument("-O", help="Catalyst with O bound", type=str)
    parser.add_argument("-CO", help="Catalyst with CO bound", type=str)
    parser.add_argument("-CN", help="Catalyst with CN bound", type=str)
    parser.add_argument("-cation", help="Ionized bare catalyst", type=str)
    parser.add_argument("-save", help="Save intermediate qdata_list pickles", type=bool, default=False)
    args = parser.parse_args() # directory names might have trailing /s

    if (os.path.isfile(args.bare+'.p') and os.path.isfile(args.cation+'.p')
        and os.path.isfile(args.O2+'_df.p') and os.path.isfile(args.OOH+'.p')
        and os.path.isfile(args.OH+'.p') and os.path.isfile(args.O+'.p')
        and os.path.isfile(args.CO+'.p')): #and os.path.isfile(args.CN+'.p')):

        cat_qdata_list = dill.load(open(args.bare + ".p", "rb"))
        cation_qdata_list = dill.load(open(args.cation + ".p", "rb"))
        O2_df = dill.load(open(args.O2 + "_df.p", "rb"))
        OOH_qdata_list = dill.load(open(args.OOH + ".p", "rb"))
        OH_qdata_list = dill.load(open(args.OH + ".p", "rb"))
        O_qdata_list = dill.load(open(args.O + ".p", "rb"))
        CO_qdata_list = dill.load(open(args.CO + ".p", "rb"))
        #CN_qdata_list = dill.load(open(args.CN + ".p", "rb"))

        for in_list in [cat_qdata_list, cation_qdata_list, OOH_qdata_list, OH_qdata_list, CO_qdata_list]:
            print(len(in_list))
    else:
        cat_qdata_list, cation_qdata_list, O2_df, OOH_qdata_list, OH_qdata_list, O_qdata_list, CO_qdata_list, CN_qdata_list = [], [], [], [], [], [], [], []
        if not cat_qdata_list == None:
            cat_qdata_list = GenLoadCat(args.bare, save=args.save)
        if not cation_qdata_list == None:
            cation_qdata_list = GenLoadCat(args.cation, save=args.save)
        if not O2_df == None:
            O2_df = GenLoadCatO2(args.O2, save=args.save)
        if not OOH_qdata_list == None:
            OOH_qdata_list = GenLoadCat(args.OOH, save=args.save)
        if not OH_qdata_list == None:
            OH_qdata_list = GenLoadCat(args.OH, save=args.save)
        if not O_qdata_list == None:
            O_qdata_list = GenLoadCat(args.O, save=args.save)
        if not CO_qdata_list == None:
            CO_qdata_list = GenLoadCat(args.CO, save=args.save)
        if not CN_qdata_list == None:
            CN_qdata_list = GenLoadCat(args.CN, save=args.save)

    match = MatchO2(O2_df, cat_qdata_list, cation_qdata_list,
        OOH_qdata_list, OH_qdata_list, O_qdata_list, CO_qdata_list) #, CN_qdata_list)
    O2_matched_df = match.catO2_df
    print(O2_matched_df)

    O2_matched_df.to_csv('catO2_matched.csv') # better to have catalyst name of first entry in file name?
