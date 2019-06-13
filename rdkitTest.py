import sys
import os
import rdkit
from rdkit import Chem




def openFiles(input_directory):    
    catalyst_dict = {}
    for file in os.listdir(input_directory):
        if file.endswith(".mol"):
            fullfilepath = input_directory + "\\" + file
            #print(file)
            catalyst_dict[file] = getNeighbors(fullfilepath)
            #print(catalyst_dict[file])
    return(catalyst_dict)
    
    
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

if __name__ == "__main__":
    input_directory = sys.argv[1]
    for file in os.listdir(input_directory):
        if file.endswith(".mol"):
            fullfilepath = input_directory + "\\" + file
            m = Chem.MolFromMolFile(fullfilepath, removeHs = False)
            if m != None:
                print(file)
                for atom in m.GetAtoms():
                    #getShortestPathToN(atom)
                    print("Index ", atom.GetIdx())
                    hasHydrogen(atom)
    #catalyst_dict = openFiles(input_directory)
    #print(catalyst_dict)
