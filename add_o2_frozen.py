# input: directory containing bare catalyst .mol and/or .xyz files
# creates .xyz files for corresponding O2-bound structures
# catalysts with one metal atom only

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math

def scale_vector(vec, new_mag):
    squared = 0.0
    new_mag = float(new_mag)
    scale = new_mag / np.linalg.norm(vec)
    return [c * scale for c in vec]

def vec3_perp(vec): # janky orthogonal vector generator; 3D
    if len(vec) != 3 or vec == [0, 0, 0]:
        print("Not a suitable vector.")
        return None
    ihat = [1, 0, 0] # arbitrarily chosen
    return scale_vector(np.cross(vec, ihat), np.linalg.norm(vec)) # retains original magnitude

indir = sys.argv[1].rstrip('/')

for infile in os.listdir(indir):
    infile_type = infile.split('.')[-1]
    if infile_type != "mol" and infile_type != "xyz":
        continue
    print(infile)

    mol3D_O2 = mol3D.mol3D()
    mol3D_O2.OBMol = mol3D_O2.getOBMol(indir + '/' + infile, infile_type, ffclean = False)
    mol3D_O2.convert2mol3D()

    if len(mol3D_O2.findMetal()) == 1:
        M_ind = mol3D_O2.findMetal()[0]
    else:
        print("%s contains %d metal atoms." %(infile, len(mol3D_O2.findMetal())))
        continue

    connection_list = mol3D_O2.getBondedAtomsSmart(M_ind, oct = False)
    #print(connection_list)

    mol, enl = structgen.ffopt('UFF', mol3D_O2, connection_list, 1, [], False, [], 200, False)

    # manually bind O2 at an angle (M-O-O = 1.8 A, 120 deg, 1.3 A)
    try:
        v1 = mol3D_O2.atoms[connection_list[0]].distancev(mol3D_O2.atoms[connection_list[1]])
        v2 = mol3D_O2.atoms[connection_list[0]].distancev(mol3D_O2.atoms[connection_list[2]])
        v = np.cross(v1, v2)
    except:
        print("Error finding plane of catalyst.")
        continue

    v_O1 = scale_vector(v, 1.8)
    v_O2 = scale_vector(v, 1.8 + 1.3 * math.sin(math.pi / 6))
    v_perp = vec3_perp(scale_vector(v, 1.3 * math.cos(math.pi / 6)))
    if v_perp == None: # shouldn't happen
        continue
    M_coords = mol3D_O2.atoms[M_ind].coords()

    mol3D_O2.addAtom(mol3D.atom3D(Sym = 'O', xyz = [M_coords[c] + v_O1[c] for c in range(len(M_coords))]))
    mol3D_O2.addAtom(mol3D.atom3D(Sym = 'O', xyz = [M_coords[c] + v_O2[c] + v_perp[c] for c in range(len(M_coords))]))

    #mol, enl = structgen.ffopt('UFF', mol3D_O2, connection_list, 1, [], False, [], 200, False)
    # i think this ffopt always fails for some reason

    mol3D_O2.writexyz(infile.split('.')[0] + "O2")
