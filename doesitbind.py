import pandas as pd
import os
import re
import sys
import getMullikan
import rdkitTest
import rdkit
from rdkit import Chem
import matplotlib.pyplot as plt
import seaborn
import getBondLengths
seaborn.set(style='ticks')


def collectData(catalyst_only_directory, O2_bound_directory, plusone_directory, molfiles_directory, xyz_directory):
    if plusone_directory != None:
        #print("here!")
        plusone_dict = {}
        
        for file in os.listdir(plusone_directory):
            #print(file)
            fullfilepath = os.path.join(plusone_directory, file)
            if file.endswith(".out"):
                #print("Processing", file)
                chelpg_part = getMullikan.getAttributes("ChElPG", fullfilepath)
            if chelpg_part == "Move on":
                continue
            cut_file = file.split('_', 1)[0]
            plusone_dict[cut_file] = chelpg_part
        
        #print(plusone_dict)
        plusone_energy_dict = {}
        
        for file in os.listdir(plusone_directory):
            #print(file)
            fullfilepath = os.path.join(plusone_directory, file)
            if file.endswith(".out"):
                #print("Processing", file)
                energy_part = getMullikan.getAttributes("Energy", fullfilepath)
            if energy_part == "Move on":
                continue
            cut_file = file.split('_', 1)[0]
            plusone_energy_dict[cut_file] = energy_part    
        #output_directory = sys.argv[3]

    if molfiles_directory != None:
        neighboring_spindensity_dict = rdkitTest.openFiles(molfiles_directory)
    
    

    catalyst_only_data = getMullikan.makeDataFrame(catalyst_only_directory)#generate catalyst only data (Mulliken spin density, ChElPG charge)
        
    unique_catalysts = catalyst_only_data.Catalyst.unique()#generate list of catalysts
    
    catalyst_dictionary = {}

    
    for catalyst in unique_catalysts:#generate dictionary of catalyst dataframes
        catalyst_dictionary[catalyst] = pd.DataFrame(columns=["Catalyst Name", "CatalystO2File", "Atom Number", "Element", "SpinDensity", "ChElPGNeutralCharge", 
        "Doesitbind", "BondLength", "BindingEnergy", "NeutralFreeEnergy","OrthoOrPara", "Meta", "FartherThanPara", "DistanceToN", 
        "AverageBondLength", "BondLengthRange", "NumberOfHydrogens", "AromaticSize", "IsInRingSize6", "IsInRingSize5", "NeighborSpinDensity", "NeighborChElPGCharge", "NeighborChargeDifference"])


    print("Individual catalyst dataframes created!")
    for index, row in catalyst_only_data.iterrows():#add data from catalyst only
        catalyst_name = row['Catalyst']
        shortened_catalyst_name = catalyst_name.split('_', 1)[0]
        dataframe = catalyst_dictionary[catalyst_name]
        appendlist = []
        appendlist.append(shortened_catalyst_name)
        appendlist.append("None")
        appendlist.append(row['Atom ID'])
        appendlist.append(row['Atoms'])
        appendlist.append(row['Mulliken Spin'])
        appendlist.append(row['ChElPG Charge'])
        appendlist.append(False)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(row['Total Free Energy'])
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)        
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        appendlist.append(0)
        dataframe.loc[index] = appendlist
    print("Catalyst data loaded!")
    
    if plusone_directory != None:
        for catalyst_name, charge in plusone_dict.items():#add data from plus one
            for catalyst in unique_catalysts:
                if catalyst_name in catalyst:
                    #print("here!")
                    dataframe = catalyst_dictionary[catalyst]
                    dataframe.insert(5, 'ChElPGPositiveCharge', charge, allow_duplicates=False)
                    
        for catalyst_name, energy in plusone_energy_dict.items():#add data from plus one
            for catalyst in unique_catalysts:
                if catalyst_name in catalyst:
                    dataframe = catalyst_dictionary[catalyst]
                    dataframe.insert(9, 'IonizedFreeEnergy', energy, allow_duplicates=False)
    
    for catalyst in unique_catalysts:
        df = catalyst_dictionary[catalyst]
        df = df.set_index("Atom Number")
        catalyst_dictionary[catalyst] = df
    
    if plusone_directory != None:
        for catalyst in unique_catalysts:
            df = catalyst_dictionary[catalyst]
            appendlist = []
            for index, row in df.iterrows():
                #print(row)
                neutral = row['ChElPGNeutralCharge']
                positive = row['ChElPGPositiveCharge']
                difference = (neutral) - float(positive)
                #print(difference)
                appendlist.append(difference)
                #print(appendlist)
            df.insert(6, 'ChargeDifference', appendlist, allow_duplicates = False)
            
        for catalyst in unique_catalysts:
            df = catalyst_dictionary[catalyst]
            appendlist = []
            for index, row in df.iterrows():
                #print(row)
                neutral = row['NeutralFreeEnergy']
                positive = row['IonizedFreeEnergy']
                difference = float(positive) - float(neutral)
                #print(difference)
                appendlist.append(difference)
                #print(appendlist)
            df.insert(10, 'IonizationEnergy', appendlist, allow_duplicates = False)    
    if molfiles_directory != None:
        for catalyst_name, dict in neighboring_spindensity_dict.items():#add neighboring spin density
            for catalyst in unique_catalysts:
                shortened_catalyst_name = catalyst.split('_', 1)[0]
                #test = catalyst.split('.', 1)[0]
                if shortened_catalyst_name in catalyst_name:
                    # print (catalyst)
                    dataframe = catalyst_dictionary[catalyst]
                    if xyz_directory != None:
                        xyz_df = getBondLengths.makeXYZ(catalyst, xyz_directory)
                        print(shortened_catalyst_name)
                        dataframe = getBondLengths.appendBondLengths(dataframe, dict, xyz_df)
                        
                    # print (dataframe.index)
                    # print (dataframe['Element'])
                    for index, row in dataframe.iterrows():
                        if dataframe.at[index, 'Element'] == 'H':
                            break
                        # print (dict)
                        listOfNeighbors = dict[index-1]
                        # print (listOfNeighbors)
                        neighborSpinDensity = 0
                        neighbor_chelpg_charge = 0
                        neighbor_charge_difference = 0
                        for neighbor in listOfNeighbors:
                            neighborSpinDensity += dataframe.at[neighbor+1, 'SpinDensity']
                            neighbor_chelpg_charge += dataframe.at[neighbor+1, 'ChElPGNeutralCharge']
                            neighbor_charge_difference += dataframe.at[neighbor+1, 'ChargeDifference']
                        dataframe.at[index, 'NeighborSpinDensity'] = neighborSpinDensity
                        dataframe.at[index, 'NeighborChElPGCharge'] = neighbor_chelpg_charge
                        dataframe.at[index, 'NeighborChargeDifference'] = neighbor_charge_difference
        for catalyst in unique_catalysts:
            dataframe = catalyst_dictionary[catalyst]
            for file in os.listdir(molfiles_directory):
                shortened_catalyst_name = catalyst.split('_', 1)[0]
                full_path = os.path.join(molfiles_directory, file)
                if shortened_catalyst_name in file:
                    m = Chem.MolFromMolFile(full_path, removeHs = False)
                    for atom in m.GetAtoms():
                        return_list = rdkitTest.getShortestPathToN(atom, file)
                        #print(return_list)
                        index1 = return_list[0]
                        #print("At index ", index1)
                        distance = return_list[1]
                        if (distance == 1):
                            dataframe.at[(index1+1), 'OrthoOrPara'] = 1
                        elif (distance == 3):
                            dataframe.at[(index1+1), 'OrthoOrPara'] = 1
                        elif distance == 2:
                            dataframe.at[(index1+1), 'Meta'] = 1
                        else:
                            dataframe.at[(index1+1), 'FartherThanPara'] = 1
                        #print("Distance is ", distance)
                        dataframe.at[(index1+1), 'DistanceToN'] = distance
                        
                        return_list_2 = rdkitTest.hasHydrogen(atom)
                        
                        index2 = return_list_2[0]
                        
                        numH = return_list_2[1]
                        dataframe.at[(index2+1), 'NumberOfHydrogens'] = numH
                        
                        return_list_3 = rdkitTest.isInRing(atom)
                        index3 = return_list_3[0]
                        ringSize = return_list_3[1]
                        if ringSize==6:
                            dataframe.at[(index3+1), 'IsInRingSize6'] = 1
                        elif ringSize==5:
                            dataframe.at[(index3+1), 'IsInRingSize5'] = 1
                            
                        return_list_4 = rdkitTest.getAromaticSize(atom)
                        index4 = return_list_4[0]
                        totalAromatic = return_list_4[1]
                        dataframe.at[(index4+1), 'AromaticSize'] = totalAromatic
                        #print (dataframe['DistanceToN'])
    
    if O2_bound_directory != None:
        O2_bound_data = pd.read_csv(O2_bound_directory)#grab O2 bound data (does it bind?)

        for index, row in O2_bound_data.iterrows():#add data from O2 binding
            catalyst_name = row['Catalyst_File_Name']
            print(catalyst_name)
            if ("\\" in catalyst_name):
                catalyst_name = catalyst_name.split("\\")[-1]
            #print(catalyst_name)
            #print(unique_catalysts)
            index = row['Active_Site']
            full_file = row['CatalystO2_File_Name']
            intended_site = re.search("-(.*)_optsp", full_file)
            intended_site = intended_site.group(1)
            intended_site = int(intended_site)
            intended_site += 1
            int_index = int(index)
            if intended_site != int_index:
                continue
            if catalyst_name in unique_catalysts:
                dataframe = catalyst_dictionary[catalyst_name]
                #print(row['BindingEnergy'])
                binding_energy = float(row['BindingEnergy'])
                bond_length = row['Cat-O2_Bond_Length']
                dataframe.at[index, 'BindingEnergy'] = binding_energy
                dataframe.at[index, 'BondLength'] = bond_length
                dataframe.at[index, 'CatalystO2File'] = full_file
                doesitbind = False
                if (bond_length<2 and binding_energy<-0.1):
                    doesitbind = True
                if dataframe.at[index, 'Doesitbind']!=True:
                    dataframe.at[index, 'Doesitbind'] = doesitbind

    print("O2 binding data loaded!")            
    alldata = pd.DataFrame()
    for catalyst in unique_catalysts:
        dataframe = catalyst_dictionary[catalyst]
        newdataframe = dataframe[(dataframe["Element"]=='C')]
        if O2_bound_directory != None:
            newdataframe = newdataframe[(newdataframe["Doesitbind"]!='Unknown')]
        alldata=alldata.append(newdataframe)
        catalyst_dictionary[catalyst] = newdataframe
    
#    for index, row in alldata.iterrows():
 #       filename = row["CatalystO2File"]
  #      short = re.search("O2-(.*)_optsp", filename)
   #     compare_part = int(short.group(1))+1
    #    int_index = int(index)
        #print (compare_part)
        #print (int_index)
     #   if compare_part != int_index:
      #      unique_id = filename
       #     alldata = alldata[(alldata["CatalystO2File"]!=filename]
            #print (alldata.shape)
    
    return (alldata)

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser("Grab data about bare catalysts and store in csv")
    parser.add_argument('-cat', help='input file/directory (bare catalyst .out)', type=str)
    parser.add_argument('-O2', help='O2 bound directory', default = None, type=str)
    parser.add_argument('-plusone', help='Plus one data for bare catalyst to grab charge difference data', default = None, type=str)
    parser.add_argument('-mol', help='Mol files for bare catalyst to grab neighbor data', default = None, type=str)
    parser.add_argument('-xyz', help='Xyz files for bare catalyst to grab bond length data', default = None, type=str)

    args = parser.parse_args()

    catalyst_only_directory = args.cat
    O2_bound_directory = args.O2
    plusone_directory = args.plusone
    molfiles_directory = args.mol
    xyz_directory = args.xyz
    
    alldata = collectData(catalyst_only_directory, O2_bound_directory, plusone_directory, molfiles_directory, xyz_directory)
    alldata.to_csv('DidItBind.csv')
    
    fg = seaborn.FacetGrid(data=alldata, hue='Doesitbind', aspect=1.61)
    fg.map(plt.scatter, 'SpinDensity', 'ChElPGNeutralCharge').add_legend()
    plt.show()
    
    #alldata.plot(kind = 'scatter', x='SpinDensity', y = 'ChElPGCharge', color = 'blue')
    #plt.show()
    