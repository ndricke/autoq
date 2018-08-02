# contains functions for randomly functionalizing atom indices
# ligand should be added to molSimplify library first
# there are currently no checks to make sure that number of requested items
    # does not exceed max possible number of items, so be careful (probably fine)

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math, random

# def find_Hs(infile): # currently unnecessary
#     # infile: .mol or .xyz
#     # returns list of indices for all Hs in molecule # might be off by one
#     infile_type = infile.split('.')[-1]
#     if not infile_type == ('mol' or 'xyz'):
#         print("Couldn't read molecule from input file.")
#         return None
#     molecule = mol3D.mol3D()
#     molecule.OBMol = molecule.getOBMol(infile, infile_type, ffclean = False)
#     molecule.convert2mol3D()
#
#     H_list = []
#     for i in range(len(molecule.atoms)):
#         if atoms[i].symbol() == 'H':
#             H_list.append(i)
#     return H_list

def find_Hs(macrocycle): # string
    # atom indices are one-indexed from .mol file
    # add custom ligands to library first, with just the metal atom removed from the .mol or .xyz
    if macrocycle == "porphyrin": # in molSimplify's library
        return [9, 10, 17, 18, 24, 25, 30, 31, 33, 34, 35, 36]
    elif macrocycle == "nan": # custom ligand (make sure these indices match the .mol you add)
        return [31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42]
    else:
        print("Not a valid ligand.")
        return None

def choose_items(list, num_items, no_repeats):
    chosen = [] # contains objects, not indices
    for n in range(num_items):
        index = random.randrange(0, len(list))
        if no_repeats: # list items should all be unique
            while list[index] in chosen:
                index = random.randrange(0, len(list))
        chosen.append(list[index])
    return chosen

def make_tempdir_outdir(macrocycle, core):
    tempdir = "/home/nricke/work/autoq/Fe-macrocycle/temp" # modify as needed
    while os.path.exists(tempdir):
        tempdir += "0"
    os.makedirs(tempdir)

    outdir = "/home/nricke/work/autoq/Fe-macrocycle/%s%s_funcs" %(macrocycle, core) # modify as needed; molSimplify expands ~ to home directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        version = 2
        outdir += str(version)
        while os.path.exists(outdir):
            version += 1
            outdir = "/home/nricke/work/autoq/Fe-macrocycle/%s%s_funcs%d" %(macrocycle, core, version)
        os.makedirs(outdir)

    return tempdir, outdir

def list_to_string(list): # no spaces
    rtn = "["
    for item in list:
        rtn += str(item) + ","
    rtn = rtn[:-1] + "]"
    return rtn

def functionalize(macrocycle, core, possible_func_indices, expected_num_funcs, num_molecules, tempdir):
    # func_library: SMILES strings or built-in ligands; modify as needed
    # generally there are issues reading in strings containing []()
    func_library = ["B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I",
                    "cyanide", "NC", "tricyanomethyl", "methylamine",
                    "dicyanamide", "nitroso",
                    "OC", "carboxyl",
                    "trifluoromethyl",
                    "thiocyanate",
                    "benzene_pi", "benzenethiol"]
    for n in range(1, num_molecules + 1):
        file_name = macrocycle + core + "-functionalized" + str(n)
        num_funcs = 0
        while (num_funcs < 1 or num_funcs > len(possible_func_indices)):
            num_funcs = np.random.poisson(expected_num_funcs)
        func_list = choose_items(func_library, num_funcs, False)
        func_indices = choose_items(possible_func_indices, num_funcs, True)

        ms_command = "molsimplify -core %s -oxstate 2 -coord 6 -geometry oct -spin 3 -lig %s -ligocc 1 -decoration %s -decoration_index %s -rundir %s -name %s" %(core, macrocycle, list_to_string(func_list), list_to_string(func_indices), tempdir, file_name)
        print(ms_command)
        os.system(ms_command)

def collect_xyz_files(current_location, destination, delete_current_location):
    # janky workaround for molSimplify's directory structure
    # copies all .xyz files in current_location to destination
    # does not check for duplicate molecules or file names
    # should this skip "badjob" files? don't expect those to happen but just in case
    for subdir, dirs, files in os.walk(current_location):
        for file in files:
            if file.split('.')[-1] == "xyz":
                os.system("scp %s %s" %(os.path.join(subdir, file), destination))
    if delete_current_location:
        os.system("ls -R %s" %current_location)
        os.system("rm -R %s" %current_location)
        os.system("rm /home/kjchen/CLIinput.inp")

def run(macrocycle, core, possible_func_indices, expected_num_funcs, num_molecules):
    tempdir, outdir = make_tempdir_outdir(macrocycle, core)
    functionalize(macrocycle, core, possible_func_indices, expected_num_funcs, num_molecules, tempdir)
    print("Functionalized .xyz files placed in " + outdir + "/")
    collect_xyz_files(tempdir, outdir, True)
    return outdir

<<<<<<< HEAD
macrocycle, core, possible_func_indices, expected_num_funcs, num_molecules = "porphyrin", "Fe", find_Hs("porphyrin"), int(sys.argv[1]), int(sys.argv[2]) # testing

tempdir, outdir = make_tempdir_outdir(macrocycle, core)

functionalize(macrocycle, core, possible_func_indices, expected_num_funcs, num_molecules)
print("Functionalized .xyz files placed in " + outdir)
collect_xyz_files(tempdir, outdir, True)
=======
if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2], find_Hs(sys.argv[1]), int(sys.argv[3]), int(sys.argv[4]))
>>>>>>> af89f49e46dd89ac93c66dae3e8b6f15afa95cc8
