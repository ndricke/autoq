import pandas as pd
import numpy as np
import sys
import os
import rdkit
from rdkit import Chem
import rdkitTest

def makeXYZ(catalyst_name, xyz_directory):#takes in catalyst_name and xyz directory, finds the relevant xyz file, returns a dataframe of all xyz coordinates for that molecule
    target = None
    catalyst_name = catalyst_name.split('.', 1)[0] #get just the name
    returndf = pd.DataFrame()
    x_list = []
    y_list = []
    z_list = []
    index_list = []
    for file in os.listdir(xyz_directory):
        if catalyst_name in file:
            target = file
    if target == None:
        print("Error occurred...")
    else:
        full_target = os.path.join(xyz_directory, target)
        with open(full_target, "r") as xyzfile:
            line = xyzfile.readline()
            line = xyzfile.readline()
            line = xyzfile.readline()
            index = 1
            while line:
                splitline = line.split()
                element = splitline[0]
                if element == 'H':
                    break
                x_list.append(float(splitline[1]))
                y_list.append(float(splitline[2]))
                z_list.append(float(splitline[3]))
                index_list.append(index)
                index+=1
                line = xyzfile.readline()

        returndf.insert(0, 'z', z_list)
        returndf.insert(0, 'y', y_list)
        returndf.insert(0, 'x', x_list)
        returndf.insert(0, 'Atom Number', index_list)
        returndf = returndf.set_index('Atom Number')
        #print(returndf)
        return(returndf)

                



def appendBondLengths(dataframe, neighbor_dict, xyz_df):#takes in dataframe to be filled in, df of xyz coordinates, and neighbor dict, and fills in dataframe with bond length data
    for index, row in xyz_df.iterrows():
        max_distance = 0
        min_distance = 100
        average_distance = 0
        listOfNeighbors = neighbor_dict[index-1]
        num_neighbors = 0
        numC = (len(xyz_df.index))
        for neighbor in listOfNeighbors:
            if neighbor >= numC:
                continue
            my_coordinates = np.array((row['x'], row['y'], row['z']))
            their_x = xyz_df.at[neighbor+1, 'x']
            their_y = xyz_df.at[neighbor+1, 'y']
            their_z = xyz_df.at[neighbor+1, 'z']
            their_coordinates = np.array((their_x, their_y, their_z))
            distance = np.linalg.norm(my_coordinates-their_coordinates)
            if distance > max_distance:
                max_distance = distance
            if distance < min_distance:
                min_distance = distance
            average_distance += distance
            num_neighbors += 1
        average_distance = average_distance/num_neighbors
        range = max_distance-min_distance
        dataframe.at[index, 'AverageBondLength'] = average_distance
        dataframe.at[index, 'BondLengthRange'] = range
    return dataframe


def getAverageBonds(dir):
    totalCCbondlengths = 0
    totalCCbonds = 0
    totalCNbondlengths = 0
    totalCNbonds = 0
    allData = pd.DataFrame()
    for file in os.listdir(dir):
        fullfilepath = os.path.join(dir, file)
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


if __name__ == "__main__":
    dir = sys.argv[1]
    makeXYZ("sf10x0_optsp_a0m2", dir)
