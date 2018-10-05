# contains class for ordered functionalization of atom indices
# ligand should be added to molSimplify library first (remove metal atom if applicable)
    # use the files on GitHub for non-porphyrin ligands
# there are currently no checks to make sure that number of requested items
    # does not exceed max possible number of items, so be careful (probably fine)

import os, sys
import numpy as np
import copy
import math, random
import itertools


import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
from autoq import decoration_manager as dm
import molSimplify.Scripts.io as msio




class FuncCatalyst(object):

    def __init__(self, catalyst, debug = False):
        # func_library: SMILES strings or built-in ligands; modify as needed
        self.catalyst = catalyst
        self.molecule, emsg = msio.lig_load(catalyst)
        self.molecule.convert2mol3D()
        self.debug = debug
        self.n = 1 #iteratively increase by 1 with each new catalyst

        #Default set of functional groups if the user doesn't input any
        self.func_library = [ "C", "N", "O", "F", "Cl", "Br", # basis set doesn't include I
                        "cyanide", "methylamine", # "NC",
                        "dicyanamide", "nitroso",
                        "OC", "carboxyl","C=O",
                        "trifluoromethyl"]
        #Default H indices, adjusting for molsimplify starting to iterate at 0. Corresponds to Ligand's index in molSimplify/Ligands
        self.Hs_dict = {"porphyrin":[9, 10, 17, 18, 24, 25, 30, 31, 33, 34, 35, 36],
                    "phencir":[31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
                    "mepyr": [18,19,24],
                    "tetry": [15,24,25,28,29,30,31,32],
                    "tetrid": [30,31,34,35,36,37,38],
                  }

    # Functionalize the catalyst at specified indices with specified functional groups
    # func_indices: contains indices to functionalized
    # func_groups: contains functional groups. Ordered to match indices
    def functionalize(self, func_indices, func_groups):
        func_lig = copy.copy(self.molecule)
        func_lig = dm.decorate_ligand(func_lig, func_groups, func_indices, debug=self.debug)
        return func_lig


    ## Functionalizes a catalyst with all possible combinations of a set of functional groups for a limited site count
    def run(self, func_num=1, func_indices=None, func_library=None):
        if func_library == None:
            func_library = self.func_library
        if func_indices == None:
            func_indices = self.Hs_dict[self.catalyst] 

        func_index_list = [item for item in itertools.combinations(func_indices, func_num)] #N Choose R for indexes to modify
        #func_group_list = [item for item in combinations(func_library, func_num)]
        func_group_list = [item for item in itertools.product(func_library, repeat=func_num)] #functional groups can repeat
        # print(func_index_list)
        # print(func_group_list)

        func_file_list = []
        for func_index_set in func_index_list:
            for func_group_set in func_group_list:

                func_molecule = self.functionalize(func_index_set, func_group_set)
                file_name = self.catalyst + "_func" + str(self.n)
                ## I modified my own custom version of molSimplify to do this, but it's not crucial for the moment
                ## Instead I'm going to just print the list of functionalizations to a file
                #func_molecule.writexyz(file_name, comment = ' '.join([self.catalyst,str(func_index_set),str(func_group_set)]))
                func_molecule.writexyz(file_name) 
                comment = ' '.join([str(self.n), self.catalyst,str(func_index_set),str(func_group_set),'\n'])
                func_file_list.append(comment)
                self.n += 1

        func_file = '_'.join([self.catalyst, str(func_num)]) + '.txt'
        with open(func_file, 'w') as f:
            for item in func_file_list:
                f.write(item)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-cat', help='Catalyst to functionalize', type=str)
    parser.add_argument('-fn', help='Number of functionalizations', type=int)
    args = parser.parse_args()

#    func_list = ['C','N','F','Cl']
#    func_list = ['C=O']

    #H_index = [24,29,32] #specific H's for tetry
    #H_index = [30, 35] #specific H's for tetrid
    H_index = None

    #func_list = None
    func_list = [ "C", "N", "O", "F", "Cl", "Br", # basis set doesn't include I
                    "cyanide", "methylamine", # "NC",
                    "OC", "C=O",
                    "trifluoromethyl"]
    func = FuncCatalyst(args.cat)
    func.run(func_num=args.fn, func_indices = H_index, func_library=func_list)
