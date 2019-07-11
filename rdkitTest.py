import pandas as pd
import numpy as np
import sys
import pandas as pd
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import collections

def getSubstructureMatch(mostCommon, molecule):
    df = pd.read_csv(mostCommon)
    list = df['SMARTS']
    #print(list)
    mol = Chem.MolFromMolFile(molecule)
    if(mol == None): return None
    subMatch = []
    for item in list:
        substructure = Chem.MolFromSmarts(item)
        #print(substructure)
        if mol.HasSubstructMatch(substructure):
            subMatch.append(1)
        else:
            subMatch.append(0)
    print(subMatch)
    return(subMatch)

def addSubstructureToCSV(csv, mostCommon, molFilesDir, numSubstructures, outdir):
    df = pd.read_csv(csv)
    numSubstructures = int(numSubstructures)
    i = 1
    while i<=numSubstructures:
        string = 'Substructure'+str(i)
        df[string] = 0
        i+=1
    for file in os.listdir(molFilesDir):
        full_file_path = os.path.join(molFilesDir, file)
        list = getSubstructureMatch(mostCommon, full_file_path)
        if list != None:
            for index, row in df.iterrows():
                if df.at[index, 'Catalyst Name'] in file:
                    i = 1
                    while i<=numSubstructures:
                        string = 'Substructure' + str(i)
                        df.at[index, string] = list[i-1]
                        print("Got 1!")
                        i += 1
    print(df)
    final_path = os.path.join(outdir, 'DidItBindWithSubstructure.csv')
    df.to_csv(final_path)
    return(df)
        

def openFiles(input_directory):    
    catalyst_dict = {}
    for file in os.listdir(input_directory):
        if file.endswith(".mol"):
            fullfilepath = os.path.join(input_directory, file)
            #print(file)
            catalyst_dict[file] = getNeighbors(fullfilepath)
            #print(catalyst_dict[file])
    return(catalyst_dict)
    
def getBondLengths(atom):
    index = atom.GetIdx()
    neighbors = atom.GetNeighbors()
    mol = atom.GetOwningMol()
    conformer = AllChem.EmbedMolecule(mol)
    conformer_obj = (mol.GetConformer(conformer))
    me_list = []
    my_id_list = []
    neighbor_list = []
    neighbor_ID_list = []
    distance_list = []
    for neighbor in neighbors:
        #print ("Begin index: ", atom.GetIdx(), "Atomic no.: ", atom.GetAtomicNum())
        #print ("Begin index: ", neighbor.GetIdx(), "Atomic no.: ", neighbor.GetAtomicNum())
        neighbor_index = neighbor.GetIdx()
        neighbor_ID = convertToElem(neighbor.GetAtomicNum())
        distance = AllChem.GetBondLength(conformer_obj, index, neighbor_index)
        me_list.append(index)
        my_id_list.append(convertToElem(atom.GetAtomicNum()))
        neighbor_list.append(neighbor_index)
        distance_list.append(distance)
        neighbor_ID_list.append(neighbor_ID)
    full_data = [me_list, my_id_list, neighbor_list, neighbor_ID_list, distance_list]
    df = pd.DataFrame(full_data)
    df = df.transpose()
    df.columns = ['Index', 'Element', 'Neighbors', 'NeighborID', 'Distance']
    return(df)
    
        
def convertToElem(num):
    if num == 6:
        return 'C'
    elif num == 7:
        return 'N'
    elif num == 1:
        return 'H'
    

    
def getNeighbors(file):
    neighbors_dict = {}
    m = Chem.MolFromMolFile(file, removeHs=False)
    if m != None:
        # print(file)
        # print(m)
        for atom in m.GetAtoms():
            index = atom.GetIdx()
            neighbors = atom.GetNeighbors()
            appendlist = []
            for atom in neighbors:
                appendlist.append(atom.GetIdx())
            neighbors_dict[index] = appendlist
        return (neighbors_dict)
        
def getNeighborsWithoutH(file):
    neighbors_dict = {}
    m = Chem.MolFromMolFile(file)
    if m != None:
        # print(file)
        # print(m)
        for atom in m.GetAtoms():
            index = atom.GetIdx()
            neighbors = atom.GetNeighbors()
            appendlist = []
            for atom in neighbors:
                appendlist.append(atom.GetIdx())
            neighbors_dict[index] = appendlist
        return (neighbors_dict)

def hasHydrogen(atom):
    #num1 = atom.GetNumExplicitHs()
    #num2 = atom.GetNumImplicitHs()
    #print("Atomic num: ", atom.GetAtomicNum())
    #print("Explicit: ", num1)
    #print("Implicit: ", num2)
    index = atom.GetIdx()
    bonds = atom.GetBonds()
    numH = 0
    for bond in bonds:
        if (bond.GetBeginAtom().GetAtomicNum() == 1) | (bond.GetEndAtom().GetAtomicNum() == 1):
            numH+=1
    #bondorder = 0
    #for bond in bonds:
     #   bondorder += bond.GetBondTypeAsDouble()
    #print("Num bonds: ", len(bonds))
    #print("Num bonds to H: ", numH)
    return_list = []
    return_list.append(index)
    return_list.append(numH)
    return return_list
    
def isInRing(atom):
    index = atom.GetIdx()
    return_list = []
    return_list.append(index)
    if atom.IsInRingSize(6):
        return_list.append(6)
        return return_list
    if atom.IsInRingSize(5):
        return_list.append(5)
        return return_list
    return_list.append(0)
    return return_list

def getAromaticSize(atom):
    return (atom.GetIdx(), len(atom.GetOwningMol().GetAromaticAtoms()))
    
def getShortestPathToN(atom, file):
    molecule = atom.GetOwningMol()
    index1 = atom.GetIdx()
    #print(index1)
    n_list = []
    for atom2 in molecule.GetAtoms():
        if atom2.GetAtomicNum() == 7:
            index2 = atom2.GetIdx()
            #print(str(atom2.GetHybridization())=="SP2")
            n_list.append(index2)
    #print (n_list)
    #print(index1)
    distance_list = []
    for value in n_list:
        if value != index1:
            #print(file)
            path = Chem.rdmolops.GetShortestPath(molecule, index1, value)
            #print("Atom is at index ", index1)
            #print("N is at index ", value)
            distance = len(path)-1
            #print(distance)
            distance_list.append(distance)
        else:
            #print("N and atom are both at index ", index1)
            return [value, 0]
    return_list = []
    return_list.append(index1)
    return_list.append(min(distance_list))
    #print(return_list)
    return return_list

def getMaxSubstructure(moldir):
    smartsList = []
    for file1 in os.listdir(moldir):
        for file2 in os.listdir(moldir):
            if file1 != file2:
                print(file1)
                print(file2)
                fullfilepath1 = os.path.join(moldir, file1)
                fullfilepath2 = os.path.join(moldir, file2)
                mol_1 = Chem.MolFromMolFile(fullfilepath1)
                mol_2 = Chem.MolFromMolFile(fullfilepath2)
                mols = [mol_1, mol_2]
                if ((mol_1 != None) & (mol_2 != None)):
                    res = rdFMCS.FindMCS(mols)
                    smarts = res.smartsString
                    smartsList.append(smarts)
    counter = collections.Counter(smartsList)
    list = counter.most_common(50)
    print(counter.most_common(50))
    df = pd.DataFrame(list)
    print(df)
    df.to_csv("D:\Kunal\Documents\MIT\\fingerprinting\mostCommon.csv")
    return df           

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Add substructure to csv")
    parser.add_argument('-csv', help='Existing csv of data from doesitbind.py', type=str)
    parser.add_argument('-mostCommon', help='csv of most common substructures in dataset', type=str)
    parser.add_argument('-mol', help='Mol files for bare catalyst to grab substructure match', default = None, type=str)
    parser.add_argument('-num', help='How many substructures to add to csv', default = None, type=str)
    parser.add_argument('-out', help='Out directory', default = None, type=str)


    args = parser.parse_args()


    addSubstructureToCSV(args.csv, args.mostCommon, args.mol, args.num, args.out)
 
    # input_directory = sys.argv[1]
    # for file in os.listdir(input_directory):
        # if file.endswith(".mol"):
            # fullfilepath = os.path.join(input_directory, file)
            # m = Chem.MolFromMolFile(fullfilepath, removeHs = False)
            # if m != None:
                # print(file)
                # for atom in m.GetAtoms():
                    # #getShortestPathToN(atom)
                    # print("Index ", atom.GetIdx())
                    # hasHydrogen(atom)
    #catalyst_dict = openFiles(input_directory)
    #print(catalyst_dict)
