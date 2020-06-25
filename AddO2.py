# input: directory containing bare catalyst .mol and/or .xyz files
# creates .xyz files for corresponding O2-bound structures
# catalysts with one metal atom only

import os, sys
import molSimplify
from molSimplify.Scripts import structgen
from molSimplify.Classes import mol3D, atom3D
import numpy as np
import math


class AddO2(object):

    def __init__(self, infile, active_sites=None):
        # also, molecule no longer explicitly optimized in this function
        self.infile = infile
        infile_type = self.infile.split('.')[-1]
        print(self.infile)
        self.OH_bond_length = 0.975

        self.molecule = mol3D.mol3D()
        self.molecule.OBMol = self.molecule.getOBMol(infile, infile_type, ffclean = False)
        self.molecule.convert2mol3D()

        if active_sites == None:
            self.active_sites = self.standardActiveSites()
        else:
            self.active_sites = active_sites

        self.bond_length_cutoffs = {'C': 1.953, 'H': 1.443, 'I': 2.543, 'Cl': 2.154,
                                   'B': 2.027, 'N': 1.834, 'F': 1.707, 'Br': 2.423,
                                   'P': 2.297, 'S': 2.198, 'O': 1.810} # X-O (angstroms)
        self.dtheta = 5 # change this if you want to
        self.binding_dict = {'O2': self.bindO2, 'O': self.bindO, 'OH': self.bindOH, 'O2H': self.bindO2H, 'CN' : self.bindCN, \
                             'CO': self.bindCO}

    def standardActiveSites(self):
        active_site_dict = {"mepyr": [14], "tetrid": [16], "tetry": [27]} # tetry[17,20]
        metal = self.molecule.findMetal()
        if len(metal) >= 1:
            self.catO_bond_length = 1.8
            self.bond_angle = 120
            self.reposition_O2 = False
            return metal
        else:
            for key in active_site_dict.keys():
                if key in self.infile:
                    self.catO_bond_length = 1.55
                    self.OO_bond_length = 1.32
                    self.bond_angle = 111
                    self.reposition_O2 = True
                    return active_site_dict[key]

    @staticmethod
    def scale_vector(vec, new_mag):
        squared = 0.0
        new_mag = float(new_mag)
        scale = new_mag / np.linalg.norm(vec)
        return [c * scale for c in vec]

    @staticmethod
    def vec3_perp(vec, theta):
        # janky orthogonal vector generator; 3D
        # theta = 0 (degrees) corresponds to i-hat direction. this is arbitrary
        if len(vec) != 3 or vec == [0, 0, 0]:
            print("Not a suitable vector.")
            return None
        v_theta = [math.cos(math.radians(theta)), math.sin(math.radians(theta)), 0]
        cross = np.cross(vec, v_theta)
        if [c for c in cross] == [0, 0, 0]: # could happen if vec and v_theta are parallel
            cross = np.cross(vec, v_theta + [0, 0, 1])
        return AddO2.scale_vector(cross, np.linalg.norm(vec)) # retains original magnitude

    def bindO2(self, site_index, molecule, v, bond_angle=111., OO_bond_length=1.32, catO_bond_length=1.55):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group O2
        """

        site_coords = molecule.atoms[site_index].coords()
        v_O1 = self.scale_vector(v, catO_bond_length)
        v_O2 = self.scale_vector(v, catO_bond_length + OO_bond_length * math.cos(math.pi - math.radians(bond_angle)))
        theta_0 = 0 # arbitrarily chosen
        v_perp = self.vec3_perp(self.scale_vector(v, self.OO_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta_0)
        if v_perp == None: # shouldn't happen
            return None
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O2[c] + v_perp[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule

    def bindO2H(self, site_index, molecule, v, bond_angle=111., OO_bond_length=1.32, catO_bond_length=1.55):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group O2H
        """

        
        site_coords = molecule.atoms[site_index].coords()
        v_O1 = self.scale_vector(v, catO_bond_length)
        v_O2 = self.scale_vector(v, catO_bond_length + OO_bond_length * math.cos(math.pi - math.radians(bond_angle)))
        v_H = self.scale_vector(v, self.OH_bond_length)
        theta_0 = 0 # arbitrarily chosen
        v_perp = self.vec3_perp(self.scale_vector(v, self.OO_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta_0)
        if v_perp == None: # shouldn't happen
            return None
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O2[c] + v_perp[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'H', xyz = [site_coords[c] + v_O2[c] + (v_perp[c]+v_H[c]) for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule

    def bindOH(self, site_index, molecule, v, bond_angle=111., OO_bond_length=1.32, catO_bond_length=1.55):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group OH
        """

        
        site_coords = molecule.atoms[site_index].coords()
        v_O1 = self.scale_vector(v, catO_bond_length)
        v_H = self.scale_vector(v, catO_bond_length + self.OH_bond_length * math.cos(math.pi - math.radians(bond_angle)))
        theta_0 = 0 # arbitrarily chosen
        v_perp = self.vec3_perp(self.scale_vector(v, self.OH_bond_length * math.sin(math.pi - math.radians(bond_angle))), theta_0)
        if v_perp == None: # shouldn't happen
            return None
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'H', xyz = [site_coords[c] + v_H[c] + v_perp[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule


    def bindCO(self, site_index, molecule, v, CO_bond_length=1.18, catC_bond_length=1.48):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group OH
        """

        
        site_coords = molecule.atoms[site_index].coords()
        v_C = self.scale_vector(v, catC_bond_length)
        v_O = self.scale_vector(v, catC_bond_length + CO_bond_length)
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'C', xyz = [site_coords[c] + v_C[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule


    def bindCN(self, site_index, molecule, v, CN_bond_length=1.16, catC_bond_length=1.49):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group OH
        """

        site_coords = molecule.atoms[site_index].coords()
        v_C = self.scale_vector(v, catC_bond_length)
        v_N = self.scale_vector(v, catC_bond_length + CN_bond_length)
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'C', xyz = [site_coords[c] + v_C[c] for c in range(len(site_coords))]))
        molecule.addAtom(mol3D.atom3D(Sym = 'N', xyz = [site_coords[c] + v_N[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule


    def bindO(self, site_index, molecule, v, catO_bond_length=1.5):
        """
        input:
            molecule: (molSimplify mol class) molecule to modify in-place with O2
            v: (np.array) coordinates to displace first O
        output:
            molecule modified with functional group O2
        """

        site_coords = molecule.atoms[site_index].coords()
        v_O1 = self.scale_vector(v, catO_bond_length)
        site_coords = molecule.atoms[site_index].coords()

        molecule.addAtom(mol3D.atom3D(Sym = 'O', xyz = [site_coords[c] + v_O1[c] for c in range(len(site_coords))]))
        #structgen.ffopt('MMFF94', molecule, [], 1, [], False, [], 200, False) # keeps failing
        return molecule

    def bindSpecies(self, site_index, species):
        """
        input:
            site_index: (list of ints) sites to bind molecule to
            species: (str) intermediate species or poison, presently limited to [O2, O2H, O, OH, CO, CN]
        output:
            None, but an xyz coordinate is saved
        Notes:
         binds O2 perpendicular to catalyst plane, then makes .xyz file
         does not modify input molecule
         bond lengths in angstroms, angle in degrees. Variables set to parameters for Fe
        """
        molecule_copy = mol3D.mol3D() #since the binding functions modify in-place, we need to make a copy here
        molecule_copy.copymol3D(self.molecule)

        connection_list = molecule_copy.getBondedAtomsSmart(site_index, oct = False)
        print("Active Site: %s       Neighbors: %s" %(str(site_index), str(connection_list))) # for debugging
        try:
            v1 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[1]])
            v2 = molecule_copy.atoms[connection_list[0]].distancev(molecule_copy.atoms[connection_list[2]])
            v = np.cross(v1, v2)
        except:
            print("Error finding plane of catalyst.")
            return None

        bindFunction = self.binding_dict[species]
        molecule_copy = bindFunction(site_index, molecule_copy, v)

        file_name_pieces = (self.infile.split('.')[0]).split('_')
        file_name = file_name_pieces[0]+species+'_'+file_name_pieces[1]
        if len(self.active_sites) > 1: #only add active site in name if more than 1
            file_name += "-" + str(site_index) # still zero-indexed
        molecule_copy.writexyz(file_name)

    def run(self, molecule):
        for site in self.active_sites:
            self.bindSpecies(site, molecule)

if __name__ == "__main__":
    infile_dir = sys.argv[1]
    species = sys.argv[2]
    if os.path.isdir(infile_dir):
        file_list = os.listdir(infile_dir)
        for infile in file_list:
            molecule_O2 = AddO2(infile_dir + '/' + infile)
            molecule_O2.run(species)
    else:
        molecule_O2 = AddO2(infile_dir)
        molecule_O2.run(species)




