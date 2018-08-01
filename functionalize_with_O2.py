# runs functionalize_catalyst.py, then add_o2_frozen.py
# creates a subdirectory of bare catalyst and another of O2-bound .xyz files

import sys, os
import functionalize_catalyst, add_o2_frozen

outdir = functionalize_catalyst.run(sys.argv[1], sys.argv[2], functionalize_catalyst.find_Hs(sys.argv[1]), int(sys.argv[3]), int(sys.argv[4]))
# macrocycle, metal/core, expected number of funcs, number of molecules generated

for file in os.listdir(outdir):
    add_o2_frozen.run(outdir + "/" + file)
os.system("mkdir %s/catfunc/ %s/catfuncO2/" %(outdir, outdir))
os.system("mv %s/*O2.xyz %s/catfuncO2/" %(outdir, outdir))
os.system("mv %s/*.xyz %s/catfunc/" %(outdir, outdir))
