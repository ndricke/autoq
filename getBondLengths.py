import pandas as pd
import numpy as np
import sys
import os
import rdkit
from rdkit import Chem
import rdkitTest

if __name__ == "__main__":
    dir = "D:\\Kunal\\Documents\\MIT\\np_catalysts\\catalystonly-molfiles"
    totalCCbondlengths = 0
    totalCCbonds = 0
    totalCNbondlengths = 0
    totalCNbonds = 0
    allData = pd.DataFrame()
    for file in os.listdir(dir):
        fullfilepath = dir + "\\" + file
        m = Chem.MolFromMolFile(fullfilepath)
        print (file)
        if m!=None:
            for atom in m.GetAtoms():
                df = rdkitTest.getBondLengths(atom)
                #print(df)
                allData = allData.append(df)
    print(allData)
    for index, row in allData.iterrows():
        if (row['Element'] == 'N') | (row['NeighborID'] == 'N'):
            totalCNbondlengths += row['Distance']
            totalCNbonds += 1
        elif (row['Element'] == 'C') & (row['NeighborID'] == 'C'):
            totalCCbondlengths += row['Distance']
            totalCCbonds += 1
    
    print("Average CN bond length: ", (totalCNbondlengths/totalCNbonds))
    print("Average CC bond length: ", (totalCCbondlengths/totalCCbonds))
                
