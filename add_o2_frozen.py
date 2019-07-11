# input: bare catalyst .mol or .xyz file/directory
# creates .xyz files for corresponding O2-bound structures
# catalysts with one metal atom only

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math
import argparse

is_metal_catalyst = True

def scale_vector(vec, new_mag):
    squared = 0.0
    new_mag = float(new_mag)
    scale = new_mag / np.linalg.norm(vec)
    return [c * scale for c in vec]

def vec3_perp(vec, theta):
    # janky orthogonal vector generator; 3D
    # theta = 0 (degrees) corresponds to i-hat direction. this is arbitrary
    if len(vec) != 3 or vec == [0, 0, 0]:
        print("Not a suitable vector.")
        return [0, 0, 0]
    v_theta = [math.cos(math.radians(theta)), math.sin(math.radians(theta)), 0]
    cross = np.cross(vec, v_theta)
    if [c for c in cross] == [0, 0, 0]: # could happen if vec and v_theta are parallel
        cross = np.cross(vec, v_theta + [0, 0, 1])
    return scale_vector(cross, np.linalg.norm(vec)) # retains original magnitude

def find_active_sites(file_name, molecule, bind_all):
    # doesn't handle catalysts with both metal and nonmetal active sites
    # use functionalize_catalyst.py first to ensure correct file names/indices
        # note that these indices are zero-indexed
    global is_metal_catalyst
    if bind_all:
        is_metal_catalyst = False
        Csp2_list = []
        molecule_copy = mol3D.mol3D()
        molecule_copy.copymol3D(molecule) # i hope this doesn't produce indexing bugs
        for i, atom in enumerate(molecule_copy.getAtoms()):
            if atom.symbol() == 'C' and len(molecule_copy.getBondedAtomsSmart(i, oct = False)) == 3:
                Csp2_list.append(i)
        print (Csp2_list)
        return Csp2_list
    if "mepyr" in file_name:
        is_metal_catalyst = False
        return [14]#[10, 12]
    elif "tetrids" in file_name:
        is_metal_catalyst = False
        return [12, 17, 24]
    elif "tetry" in file_name:
        is_metal_catalyst = False
        return [21, 22, 23]
    if len(molecule.findMetal()) > 0:
        is_metal_catalyst = True # i think this is redundant
        return molecule.findMetal()
    else:
        print("Could not identify the active site(s) in " + file_name)
        return None

def get_normal_vec(molecule_copy, site_index):
    connection_list = molecule_copy.getBondedAtomsSmart(site_index, oct = False)
    #print("Active Site: %s       Neighbors: %s" %(str(site_index), str(connection_list))) # for debugging
    try:
        v1 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[1]])
        v2 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[2]])
        v = np.cross(v1, v2)
        return v
    except:
        print("Error finding plane of catalyst.")
        return None # the script doesn't always stop when this happens so be careful

# bind_OOH() doesn't check for collisions but if this turns
    # out to be an actual issue then I'll add that in

def bind_OOH(infile, molecule, site_index, catO_bond_length, catOO_bond_angle, OO_bond_length, OH_bond_length):
    molecule_copy = mol3D.mol3D()
    molecule_copy.copymol3D(molecule)

    v = get_normal_vec(molecule_copy, site_index)
    v_O1 = scale_vector(v, catO_bond_length)
    v_O2 = scale_vector(v, catO_bond_length + OO_bond_length * math.cos(math.pi - math.radians(catOO_bond_angle)))
    theta_0 = 0
    v_perp = vec3_perp(scale_vector(v, OO_bond_length * math.sin(math.pi - math.radians(catOO_bond_angle))), theta_0)
    if v_perp == None: # shouldn't happen
        return None
    v_H = scale_vector(v, OH_bond_length)

    site_coords = molecule_copy.atoms[site_index].coords()

    molecule_copy.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
    molecule_copy.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O2[c] + v_perp[c] for c in range(len(site_coords))]))
    molecule_copy.addAtom(mol3D.atom3D(Sym = 'H', xyz = [site_coords[c] + v_O2[c] + v_perp[c] + v_H[c] for c in range(len(site_coords))]))

    file_name = infile.split('.')[0] + "OOH"
    if not is_metal_catalyst: # metal catalysts assumed to have one active site
        file_name += "-" + str(site_index) # still zero-indexed
    molecule_copy.writexyz(file_name)

def bind_O(infile, molecule, site_index, catO_bond_length):
    molecule_copy = mol3D.mol3D()
    molecule_copy.copymol3D(molecule)

    v = get_normal_vec(molecule_copy, site_index)
    v_O = scale_vector(v, catO_bond_length)

    site_coords = molecule_copy.atoms[site_index].coords()

    molecule_copy.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O[c] for c in range(len(site_coords))]))

    file_name = infile.split('.')[0] + "O"
    if not is_metal_catalyst: # metal catalysts assumed to have one active site
        file_name += "-" + str(site_index) # still zero-indexed
    molecule_copy.writexyz(file_name)

def bind_XY(infile, molecule, site_index, X, Y, catX_bond_length, bond_angle, XY_bond_length, reposition_XY, fn_suffix):
    # XY is any diatomic species, with X closer to the active site
    # positions XY perpendicular to catalyst plane, then makes .xyz file
        # does not modify input molecule
    # bond lengths in angstroms, angle in degrees
    molecule_copy = mol3D.mol3D()
    molecule_copy.copymol3D(molecule)

    v = get_normal_vec(molecule_copy, site_index)
    v_X = scale_vector(v, catX_bond_length)
    v_Y = scale_vector(v, catX_bond_length + XY_bond_length * math.cos(math.pi - math.radians(bond_angle)))
    theta_0 = 0 # arbitrarily chosen
    v_perp = vec3_perp(scale_vector(v, XY_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta_0)
    if v_perp == None: # shouldn't happen
        return None

    site_coords = molecule_copy.atoms[site_index].coords()

    molecule_copy.addAtom(mol3D.atom3D(Sym = X, xyz = [site_coords[c] + v_X[c] for c in range(len(site_coords))]))
    molecule_copy.addAtom(mol3D.atom3D(Sym = Y, xyz = [site_coords[c] + v_Y[c] + v_perp[c] for c in range(len(site_coords))]))
    #structgen.ffopt('MMFF94', molecule_copy, [], 1, [], False, [], 200, False) # keeps failing

    successful_bind = True
    if reposition_XY: # workaround for ffopt() failing to set up; nonmetal catalysts
        bond_length_cutoffs = {'C': 1.953, 'H': 1.443, 'I': 2.543, 'Cl': 2.154,
                               'B': 2.027, 'N': 1.834, 'F': 1.707, 'Br': 2.423,
                               'P': 2.297, 'S': 2.198, 'O': 1.810} # bound to O (angstroms); X, Y may not be O but oh well
        dtheta = 5 # change this if you want to
        theta = theta_0 + dtheta
        position_OK, original_direction = False, True
        while not position_OK:
            position_OK = True
            for atom in molecule_copy.getAtoms():
                if atom == molecule_copy.getAtoms()[-1] or atom == molecule_copy.getAtoms()[-2]:
                    continue
                if atom.symbol() in bond_length_cutoffs:
                    cutoff = bond_length_cutoffs[atom.symbol()]
                else: # won't happen if bond_length_cutoffs is kept up to date
                    cutoff = bond_length_cutoffs[max(bond_length_cutoffs, key = lambda key: bond_length_cutoffs[key])]
                    print("Using %f A as the bond length cutoff for %s" %(cutoff, atom.symbol()))
                if molecule_copy.getAtoms()[-1].distance(atom) < cutoff:
                    position_OK = False
                    break
            if not position_OK:
                if theta >= 360:
                    if not original_direction:
                        successful_bind = False
                        print("%s may have geometry issues after O2-binding." %infile)
                        break
                    original_direction = False
                    v_X = [-c for c in v_X]
                    v_Y = [-c for c in v_Y]
                    theta = 0
                v_perp = vec3_perp(scale_vector(v, XY_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta)
                if v_perp == None: # shouldn't happen
                    return None
                molecule_copy.getAtoms()[-2].setcoords([site_coords[c] + v_X[c] for c in range(len(site_coords))]) # X
                molecule_copy.getAtoms()[-1].setcoords([site_coords[c] + v_Y[c] + v_perp[c] for c in range(len(site_coords))]) # Y
                theta += dtheta
                print("reposition", infile, site_index, original_direction)

    file_name = infile.split('.')[0] + fn_suffix
    if not is_metal_catalyst: # metal catalysts assumed to have one active site
        file_name += "-" + str(site_index) # still zero-indexed
    molecule_copy.writexyz(file_name)

    if not successful_bind:
        error_outloc = os.path.dirname(infile)
        if error_outloc == "":
            error_outloc = os.getcwd()
        error_fn = "bind_error_log.txt"
        with open(error_outloc + '/' + error_fn, 'a') as error_file:
            error_file.write("%s|%d\n" %(file_name, site_index)) # might create duplicate entries

def run(infile, add_O2, add_O2_reactant, O2r_dist, add_OOH_O_OH, add_CO_CN, bind_all):
    # changed from indir to infile to prevent segfaults
        # also, molecule no longer explicitly optimized in this function
    infile_type = infile.split('.')[-1]
    if infile_type != "mol" and infile_type != "xyz":
        return None
    print(infile)

    mol3D_O2 = mol3D.mol3D()
    mol3D_O2.OBMol = mol3D_O2.getOBMol(infile, infile_type, ffclean = False)
    mol3D_O2.convert2mol3D()

    active_sites = find_active_sites(infile, mol3D_O2, bind_all)
    if active_sites == None:
        return None
    if is_metal_catalyst:
        for site in active_sites: # expect len(active_sites) to be 1; don't set bind_all == True
            if add_O2:
                bind_XY(infile, mol3D_O2, site, "O", "O", 1.8, 120, 1.3, False, "O2")
            if add_O2_reactant:
                bind_XY(infile, mol3D_O2, site, "O", "O", O2r_dist, 120, 1.21, False, "-O2-%2.3f" %O2r_dist)
            if add_OOH_O_OH:
                bind_OOH(infile, mol3D_O2, site, 1.762, 114.1, 1.457, 0.978)
                bind_O(infile, mol3D_O2, site, 1.619)
                bind_XY(infile, mol3D_O2, site, "O", "H", 1.852, 119.7, 0.973, False, "OH")
            if add_CO_CN:
                bind_XY(infile, mol3D_O2, site, "C", "O", 1.701, 180, 1.161, False, "CO")
                bind_XY(infile, mol3D_O2, site, "C", "N", 1.847, 180, 1.179, False, "CN")
    else: # all functions already work for nonmetal catalysts-- just need the correct bond lengths
        for site in active_sites:
            if add_O2:
                bind_XY(infile, mol3D_O2, site, "O", "O", 1.55, 111, 1.3, True, "O2")
            if add_O2_reactant:
                bind_XY(infile, mol3D_O2, site, "O", "O", O2r_dist, 120, 1.21, False, "-O2-%2.3f" %O2r_dist)
            if add_OOH_O_OH:
                bind_OOH(infile, mol3D_O2, site, 1.445, 107.3, 1.462, 0.978)
                bind_O(infile, mol3D_O2, site, 1.364)
                bind_XY(infile, mol3D_O2, site, "O", "H", 1.433, 107.0, 0.975, False, "OH")
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='input file/directory (bare catalyst .mol or .xyz)', type=str)
    parser.add_argument('-O2', help='binds O2 to catalyst', type=bool,default=False) # parser parses first argument incorrectly sometimes?
    parser.add_argument('-O2r', help='adds unbound O2 near catalyst', type=bool,default=False)
    parser.add_argument('-intermediate', help='binds OOH, O, OH to catalyst', type=bool,default=False)
    parser.add_argument('-poison', help='binds CO, CN to catalyst', type=bool,default=False)
    parser.add_argument('-O2rdist', help='distance (angstroms) between active site and unbound O2', type=float,default=3) # default should be enough to ensure that O2 not bound to catalyst
    parser.add_argument('-bindAll', help='identifies all/only sp^2 carbons as active sites (assumes nonmetal catalyst)', type=bool,default=False)
    args = parser.parse_args()

    if os.path.isfile(args.f):
        run(args.f, args.O2, args.O2r, args.O2rdist, args.intermediate, args.poison, args.bindAll)
    elif os.path.isdir(args.f):
        files = []
        for file in os.listdir(args.f):
            if file.split('.')[-1] == "xyz" or file.split('.') == "mol":
                files.append(file)
        for file in files: # fixes issues caused by adding new files to a directory you're iterating over
            run(args.f.rstrip('/') + '/' + file, args.O2, args.O2r, args.O2rdist, args.intermediate, args.poison, args.bindAll)
    else:
        print("Invalid input for -f")