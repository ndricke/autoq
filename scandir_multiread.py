import sys
import os
import argparse
import pandas as pd
import numpy as np
import re
from ast import literal_eval as make_tuple

import scandir_analysis

"""
This script is for reading a lot of scandir output json or csv's
This came up for the autoq ngcc project where all the files were in cat-out, cat-CO-out, cat-O2H-out... directories
The reason I couldn't combine all of the output files into one directory is that the naming structure for functionalizing
restarted from 1 for each directory, so there were repeats

This script needs to read each file, and use the various base parts of names to re-map the numbering for all of the catalysts
(I don't technically need the mapping of functional groups to data points for the pair-plot or the volcano plot, but I do
need to be able to link all the various intermediates together, which means I can't have name collisions)

I could just rename everything with a script on it's own, which also operates on the functionalization mappings. This is also going to be necessary for the macrocycles.



"""

def load_funclist(fname):
    funcnum_list, cat_list, loc_list, func_list = [], [], [], []
    with open(fname, "r") as f:
        for line in f:
            spline = line.split()
            funcnum_list.append(int(spline[0]))
            cat_list.append(spline[1])
            loc_list.append(make_tuple(spline[2]))
            func_list.append(make_tuple(spline[3]))
    return pd.DataFrame({"funcnum":funcnum_list, "catalyst":cat_list, "loc":loc_list, "func":func_list})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir", help="Input data directory containing scandir jsons", default="scandir_jsons")
    parser.add_argument("-f", "--funclist", help="Input directory containing functionalization lists", default="catfunc_list")
    parser.add_argument("-o", "--outfile", help="Output file with concatenated data", default="catdata_all.json")
    args = parser.parse_args()

    print("Parsing functionalization map files")
    func_map = {}
    for funclist_filename in os.listdir(args.funclist):
        print(funclist_filename)
        funclist_filepath = args.funclist + '/' + funclist_filename
        funclist = funclist_filename.split('.')[0]
        func_map[funclist] = load_funclist(funclist_filepath)
    print()

    print("Reading and analyzing scandir jsons")
    df_list = []
    for infile in os.listdir(args.datadir):
        fpath = args.datadir + '/'+infile
        f_prefix = infile.split(".")[0]
        df = pd.read_json(fpath)
        df_aug_min = scandir_analysis.parse_autoq_catalysts(df)
        df_aug_min["data_dir"] = f_prefix
        df_list.append(df_aug_min)
    df_all = pd.concat(df_list)
    df_all.reset_index(inplace=True)
    df_all.to_json(args.outfile)


