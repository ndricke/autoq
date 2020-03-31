import pandas as pd
import os
import re
import sys
import matplotlib.pyplot as plt
import seaborn


def makeEnergyDeltaDf(df):#dataframe comes from energies by active site cleaned
    deltaO2list = []
    deltaOOHlist = []
    deltaOlist = []
    deltaOHlist = []
    
    for index, row in df.iterrows():
        deltaO2list.append((float(row['O2Energy']))-(float(row['catEnergy'])))
        deltaOOHlist.append((float(row['OOHEnergy']))-(float(row['O2Energy'])))
        deltaOlist.append((float(row['OEnergy']))-(float(row['OOHEnergy'])))
        deltaOHlist.append((float(row['OHEnergy']))-(float(row['OEnergy'])))
    
    dict = {
    'deltaO2': deltaO2list,
    'deltaOOH': deltaOOHlist,
    'deltaO': deltaOlist,
    'deltaOH': deltaOHlist
    }
    df = pd.DataFrame(dict)
    print(df)
    return(df)