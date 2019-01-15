# runs functionalize_catalyst.py, then add_o2_frozen.py
# creates a subdirectory of bare catalyst and another of O2-bound .xyz files
# see functionalize_catalyst.find_Hs() for list of accepted catalysts

import sys, os
import functionalize_catalyst, add_o2_frozen
import argparse

parser = argparse.ArgumentParser()
# for functionalize_catalyst
parser.add_argument('-cat', help='Catalyst name (porphyrin, nan, mepyrid, tetrids, or tetry)', type=str)
parser.add_argument('-core', help='Metal atom', type=str,default="Fe")
parser.add_argument('-numFunc', help='Expected number of functionalizations per molecule', type=int)
parser.add_argument('-numMol', help='Number of molecules generated', type=int)

# for add_o2_frozen
parser.add_argument('-O2', help='binds O2 to catalyst', type=bool,default=True)
parser.add_argument('-O2r', help='adds unbound O2 near catalyst', type=bool,default=False)
parser.add_argument('-intermediate', help='binds OOH, O, OH to catalyst', type=bool,default=False)
parser.add_argument('-poison', help='binds CO, CN to catalyst', type=bool,default=False)
args = parser.parse_args()

outdir = functionalize_catalyst.run(args.cat, args.core, functionalize_catalyst.find_Hs(args.cat), args.numFunc, args.numMol)
# for nonmetal catalysts, you can put anything for the core (functionalize_catalyst.py ignores it)

for file in os.listdir(outdir):
    add_o2_frozen.run(outdir + "/" + file, args.O2, args.O2r, args.intermediate, args.poison)

# won't work properly if catalyst name contains "O2", etc.
if args.O2r:
    os.system("mkdir %s/catfunc-O2/" %outdir)
    os.system("mv %s/*-O2*.xyz %s/catfunc-O2/" %(outdir, outdir))
if args.O2:
    os.system("mkdir %s/catfuncO2/" %outdir)
    os.system("mv %s/*O2*.xyz %s/catfuncO2/" %(outdir, outdir))
if args.poison:
    os.system("mkdir %s/catfuncCO/ %s/catfuncCN/" %(outdir, outdir))
    os.system("mv %s/*CO*.xyz %s/catfuncCO/" %(outdir, outdir))
    os.system("mv %s/*CN*.xyz %s/catfuncCN/" %(outdir, outdir))
if args.intermediate:
    os.system("mkdir %s/catfuncOOH/ %s/catfuncO/ %s/catfuncOH/" %(outdir, outdir, outdir))
    os.system("mv %s/*OOH*.xyz %s/catfuncOOH/" %(outdir, outdir))
    os.system("mv %s/*OH*.xyz %s/catfuncOH/" %(outdir, outdir))
    os.system("mv %s/*O*.xyz %s/catfuncO/" %(outdir, outdir))
os.system("mkdir %s/catfunc/" %outdir)
os.system("mv %s/*.xyz %s/catfunc/" %(outdir, outdir))
