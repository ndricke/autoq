# runs functionalize_catalyst.py, then add_o2_frozen.py
# creates a subdirectory of bare catalyst and another of O2-bound .xyz files
# see functionalize_catalyst.find_Hs() for list of accepted catalysts

import sys, os
import functionalize_catalyst, add_o2_frozen

outdir = functionalize_catalyst.run(sys.argv[1], sys.argv[2], functionalize_catalyst.find_Hs(sys.argv[1]), int(sys.argv[3]), int(sys.argv[4]))
# catalyst, metal/core, expected number of funcs, number of molecules generated
# for nonmetal catalysts, you can put anything for the core (functionalize_catalyst.py ignores it)

for file in os.listdir(outdir):
    add_o2_frozen.run(outdir + "/" + file, True, True)
    # last two parameters control which species are bound to catalyst

os.system("mkdir %s/catfunc/ %s/catfuncO2/ %s/catfuncOOH/ %s/catfuncO/ %s/catfuncOH/" %(outdir, outdir, outdir, outdir, outdir))
os.system("mv %s/*O2*.xyz %s/catfuncO2/" %(outdir, outdir))
os.system("mv %s/*OOH*.xyz %s/catfuncOOH/" %(outdir, outdir))
os.system("mv %s/*OH*.xyz %s/catfuncOH/" %(outdir, outdir))
os.system("mv %s/*O*.xyz %s/catfuncO/" %(outdir, outdir))
# won't work properly if catalyst name contains "O2", etc.
os.system("mv %s/*.xyz %s/catfunc/" %(outdir, outdir))
