import pandas as pd
import numpy as np
import os
import sys
import re

class CatalystActiveSite:
    def __init__(self, catalyst, activesite):
        self.catalyst = catalyst
        self.activesite = activesite
        self.catEnergy = 0
        self.OEnergy = 0
        self.O2Energy = 0
        self.OHEnergy = 0
        self.OOHEnergy = 0

def reorganizeIntoCatalysts(df):#csv comes from organizeEnergies
    objectlist = []
    tuplelist = []
    for index, row in df.iterrows():
        catalyst = row['Catalyst']
        activesite = row['ActiveSite']
        if activesite != 'Cat':
            newtuple = (catalyst,activesite)
            tuplelist.append(newtuple)
    y = np.unique(tuplelist, axis=0)
    z = [] 
    for i in y:
        z.append(i)
    for tuple in z:
        x = CatalystActiveSite(tuple[0], tuple[1])
        print(x.catalyst)
        print(x.activesite)
        objectlist.append(x)
    for object in objectlist:
        for index, row in df.iterrows():
            if row['Catalyst'] == object.catalyst:
                if row['Intermediate'] == 'Cat':
                    object.catEnergy = row['Energy']
                if row['ActiveSite'] == object.activesite:
                    if row['Intermediate'] == 'O2':
                        object.O2Energy = row['Energy']
                    elif row['Intermediate'] == 'OOH':
                        object.OOHEnergy = row['Energy']
                    elif row['Intermediate'] == 'O':
                        object.OEnergy = row['Energy']
                    elif row['Intermediate'] == 'OH':
                        object.OHEnergy = row['Energy']

    df = pd.DataFrame([vars(object) for object in objectlist])
    print(df)
    return(df)

def collectEnergies(filename):
    with open(filename, "r") as file:
        total_free_energy = 0
        line = file.readline()
        keyphrase = "Total Free Energy"
        deathphrase = "Q-Chem fatal error"
        while line:
            if deathphrase in line:
                #print("This file has a fatal error! Moving on...")
                return(0)
            if keyphrase in line:
                splitline = line.split()
                total_free_energy = float(splitline[9])
            line = file.readline()
        return total_free_energy
        
def organizeEnergies(indir):
    df = pd.DataFrame()
    filelist = []
    energylist = []
    intermediatelist = []
    catalystlist = []
    activesitelist = []
    for file in os.listdir(indir):
        full_path = os.path.join(indir, file)
        filelist.append(file)
        energylist.append(collectEnergies(full_path))
        if 'O' in file:
            catalyst = file.split('O')[0].split('f')[1]
        else:
            catalyst = file.split('_')[0].split('f')[1]
        catalystlist.append(catalyst)
        try:
            intermediate = file.split('-')[0].split('x')[1][1:]
            if len(intermediate) > 3:
                intermediate = 'Cat'
        except: 
            intermediate = 'Cat'   
        intermediatelist.append(intermediate)
        try:
            activesite = file.split('-')[1].split('_')[0]
        except:
            activesite = 'Cat'
        activesitelist.append(activesite)
        print(file)
    df['Filename'] = filelist
    df['Energy'] = energylist
    df['Intermediate'] = intermediatelist
    df['Catalyst'] = catalystlist
    df['ActiveSite'] = activesitelist
    return(df)

if __name__ == "__main__":
    df = organizeEnergies(sys.argv[1])
    outdir = sys.argv[2]
    csv_name = 'intermediateEnergies.csv'
    outpath = os.path.join(outdir, csv_name)
    print(df)
    df.to_csv(outpath)
    print("To csv!")
    newdf = reorganizeIntoCatalysts(df)
    reorganized_csv_name = 'energiesByActiveSite.csv'
    outpath = os.path.join(outdir, reorganized_csv_name)
    newdf.to_csv(outpath)
    print("To csv!")