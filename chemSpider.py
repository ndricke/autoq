import pandas as pd
import os
import sys
import re
import numpy as np 
  
# function to get unique values 
def makeUnique(list1): 
    print("here!")
    x = np.array(list1) 
    return(np.unique(x))

def parseSpider(filename):
    smileslist = []
    targetString = "SMILES"
    with open(filename, "r") as file:
        line = file.readline()
        while line:
            if targetString in line:
                chemformula = line.split('SMILES')[1].split('.')
                #print (chemformula)
                #print(chemformula)
                for ion in chemformula:
                    ion = ' '.join(ion.split())
                    if "n+" in ion or "N+" in ion:
                        #print("n+")
                        #print(ion)
                        smileslist.append(ion)
                        #print("I appended ", ion)
            line = file.readline()
    uniquelist = makeUnique(smileslist)
    return(uniquelist)
    
def makeSmiles(smileslist):
    counter = 39
    for smiles in smileslist:
        filename = "D:\Kunal\Documents\MIT\chemspider"
        charge = smiles.count("n+") + smiles.count("N+")
        if charge == 1:
            filename = filename + "\smiles_neutral\sf" + str(counter) + "x0_optsp_" + "a0m2.smi"
        elif charge == 2:
            filename = filename + "\smiles_c1\sf" + str(counter) + "x0_optsp_" + "c1m2.smi"
        elif charge == 3:
            filename = filename + "\smiles_c2\sf" + str(counter) + "x0_optsp_" + "c2m2.smi"
        else:
            print("Found an odd charge")
            continue
        f = open(filename, "w+")
        f.write(smiles)
        f.close()
        counter += 1

if __name__ == "__main__":
    returnlist = parseSpider(sys.argv[1])
    print(returnlist)
    makeSmiles(returnlist)