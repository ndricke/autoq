import os
import re
import sys
import shutil

CuDir = r"D:\Kunal\Documents\MIT\metal_centered_catalysis\CuTMPA\out_xyz"

samplefile = "D:\Kunal\Documents\MIT\metal_centered_catalysis\CuTMPA\out_xyz\2_optsp_c1m1.xyz"
def getExtraCharge(file):
    if '5' in file:
        return 1
    elif '6' in file:
        return 1
    elif '9' in file:
        return 1
    elif '8' in file:
        return 2
    else:
        return 0

def getExtraElectrons(file):
    if '5' in file:
        return 1
    elif '6' in file:
        return 1
    elif '9' in file:
        return 1
    elif '8' in file:
        return 2
    else:
        return 0

def IsBimetallic(file):
    if '3' in file:
        return True
    elif '6' in file:
        return True
    elif '8' in file:
        return True
    else:
        return False
        
def getMetalElectronCount(metal):
    metaldict = {
    'Mn' : 25, 
    'Fe' : 26,
    'Co' : 27,
    'Ni' : 28,
    'Cu' : 29
    }
    return metaldict[metal]    

def getMetalChargeStates(metal):
    chargedict = {
    'Mn' : (1,2,3,4), 
    'Fe' : (1,2,3),
    'Co' : (1,2,3),
    'Ni' : (1,2),
    'Cu' : (1,2)    
    }
    return chargedict[metal]

def swapMetal(srcfile, destdir, destfilename, metal):
    final_path = os.path.join(destdir, destfilename)
    shutil.copyfile(srcfile, final_path)
    with open(final_path, 'r') as file:
        data = file.readlines()
    
    for i in range(len(data)):
        if 'Cu' in data[i]:
            data[i] = data[i].replace('Cu', metal)

    with open(final_path, 'w') as file:
        file.writelines(data)
        
def createNewMetalDirectory(metal):
    dir = r"D:\Kunal\Documents\MIT\metal_centered_catalysis"
    for file in os.listdir(CuDir):
        if file.endswith(".xyz"):
            srcfile = os.path.join(CuDir, file)
            for oxidationstate in getMetalChargeStates(metal):
                destfilename = file[0] + "_optsp_"
                charge = oxidationstate + getExtraCharge(file)
                eCount = getMetalElectronCount(metal)-oxidationstate
                eCount += getExtraElectrons(file)
                if IsBimetallic(file):
                    eCount = (eCount*2)
                if (eCount % 2) == 0:
                    mult = 1
                if (eCount % 2) == 1:
                    mult = 2
                if mult == 1:
                    finalfilename = destfilename + "c" + str(charge) + "m" + "1.xyz"
                    folder = metal + "TMPA"
                    destdir = os.path.join(dir, folder)
                    swapMetal(srcfile, destdir, finalfilename, metal)
                    print("Created ", finalfilename)
                    finalfilename = destfilename + "c" + str(charge) + "m" + "3.xyz"
                    swapMetal(srcfile, destdir, finalfilename, metal)
                    print("Created ", finalfilename)
                elif mult == 2:
                    destfilename = destfilename + "c" + str(charge) + "m" + "2.xyz"
                    folder = metal + "TMPA"
                    destdir = os.path.join(dir, folder)
                    swapMetal(srcfile, destdir, destfilename, metal)
                    print("Created ", destfilename)

if __name__ == "__main__":
    createNewMetalDirectory("Fe")
    createNewMetalDirectory("Co")
    createNewMetalDirectory("Ni")
    createNewMetalDirectory("Mn")

