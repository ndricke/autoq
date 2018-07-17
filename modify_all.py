
import os
import pandas as pd
from rdkit import Chem
from Modifications2 import random_on_target
from Modifications2 import random_from_list
from Modifications2 import weighted
from shutil import copy
import argparse
from random import choice

parser = argparse.ArgumentParser()
parser.add_argument("--pardir", type=str, help="the parent directory from which to select catalysts to be modified")
parser.add_argument("--moddir", type=str, help="the name of the destination directory for the modified catalysts")
parser.add_argument("--modnum", type=int, help="number of modifications to be applied to each catalyst")
args = parser.parse_args()


mol_dir=args.pardir
mod_dir = args.moddir
modif_no=args.modnum

print(mol_dir)
path = mol_dir+'/'

cutoff=7

for c_name in os.listdir(mol_dir):
    print('\n'+c_name+'\n\n')
    
    c_mol_file=path + c_name 
    m=Chem.MolFromMolFile(c_mol_file, removeHs=False)
    mod_list=[Chem.MolToSmiles(m)]     
    at_list=list(range(cutoff+1,len(m.GetAtoms())))
    n=general_modify(m,cutoff)
    
    for j in range(modif_no):
        print('working on modif '+str(j+1))
        
        if n!=None:
        
            while Chem.MolToSmiles(n) in mod_list:
                n=random_on_target(m,cutoff)
                
            mod_list.append(Chem.MolToSmiles(n))
            
            xyz_file=mod_dir+'/'+c_name+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.xyz'
            mol_file=mod_dir+'/'+c_name+'_m'+(2-len(str(j+1)))*'0'+str(j+1)+'.mol'
            

            #write mol file for bare modified catalyst
            with open(mol_file,'w') as g:
                g.write(Chem.MolToMolBlock(n))
                
            #mol to xyz
            os.system("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
#                print("babel -imol %s -oxyz %s" % (mol_file,xyz_file))
            
            ##run add_o2_frozen which creates .xyz and .mol files for the bare and o2 bound structures
            #os.system('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
            #print('python add_o2_frozen.py '+mod_c_dir+' '+str(cutoff))
                

        else:
            print('cannot modify '+c_name+' anymore')

