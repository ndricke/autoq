import sys
import os
import copy



indir = sys.argv[1]
outdir = sys.argv[2]

for infile in os.listdir(indir):
    print("Starting %s" % infile)
    H_lines = []
    with open(indir+"/"+infile, 'r') as f:
        fl = list(f)
        for i, line in enumerate(fl):
            if line.split()[0] == 'H':
                H_lines.append(i)

    print(H_lines)
    print(len(fl))
    reorder_H = [fl[i] for i in range(len(fl)) if i not in H_lines]

    for H_index in H_lines:
        reorder_H.append(fl[H_index])

    outfile = outdir+"/"+infile
    with open(outfile, 'w') as outf:
        for item in reorder_H:
            outf.write(item)

    print("Done with %s" % infile)
 





