import os
import sys
import shutil
import pandas as pd

def getSitesFromCSV(csv, indir, outdir):
    df = pd.read_csv(csv)
    df = df[df["Doesitbind"]==True]
    list =  df["CatalystO2File"].tolist()
    #print(list)
    src_files = os.listdir(indir)
    for file in src_files:
        full_file_path = os.path.join(indir, file)
        #print(full_file_path)
        if full_file_path in list:
            print(full_file_path)
            shutil.copy(full_file_path, outdir)
            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Get binding sites")
    parser.add_argument('-csv', help = 'doesitbind output csv', type=str)
    parser.add_argument('-indir', help='input directory', type=str)
    parser.add_argument('-outdir', help='output directory', type=str)
    
    args = parser.parse_args()

    getSitesFromCSV(args.csv, args.indir, args.outdir)