# input: reactant (unbound catalyst) .xyz file/directory and product (O2-bound) .xyz file/directory
# matches reactant to product file and writes an fsm .in file

import argparse
import os

# same as O2_binding.catalyst_name()
def catalyst_name(str): # returns name of bare catalyst given catalystO2 file name
    str_bn = os.path.basename(str)
    return str_bn[0:str_bn[0:str_bn.find("_")].rfind("O2")] # basename

def find_reactant_file(product_file, reactant_dir):
    for file in os.listdir(reactant_dir):
        if file[0:len(catalyst_name(product_file))] == catalyst_name(product_file) and file.split('.')[-1] == "xyz":
            return file # basename
    return None

def write_fsm_file(reactant_file, product_file, c, m, basis, method, node): # reactant_file and product_file are paths
    # writes fsm file to same directory as product_file to reduce risk of overwriting
    output_location = os.path.dirname(product_file)
    if output_location == "":
        output_location = os.getcwd()

    charge_type = "a" # does the same thing as charge_tran in qcIn.py
    if c > 0:
        charge_type = "c"
    product_file_bn = os.path.basename(product_file)
    fsm_fn = "%s_%s_%s%d%s%d.in" %(product_file_bn[0:product_file_bn.rfind('.')].split('_')[0], 'fsm', charge_type, c, 'm', m)

    fsm_file = open(output_location + '/' + fsm_fn, 'w')

    fsm_file.write("$molecule\n")
    fsm_file.write("%d %d\n" %(c, m))

    with open(reactant_file) as open_reactant_file:
        for line_num, line in enumerate(open_reactant_file, start = 1):
            if line_num >= 3 and line.strip() != "":
                fsm_file.write(line.strip() + "\n")

    fsm_file.write("****\n")

    with open(product_file) as open_product_file:
        for line_num, line in enumerate(open_product_file, start = 1):
            if line_num >= 3 and line.strip() != "":
                fsm_file.write(line.strip() + "\n")

    fsm_file.write("$end\n\n")
    fsm_file.write("$rem\n")
    fsm_file.write("JOBTYPE                       fsm\n")
    fsm_file.write("FSM_NGRAD                     4\n")
    fsm_file.write("FSM_NNODE                     %d\n" %node)
    fsm_file.write("FSM_MODE                      2\n")
    fsm_file.write("FSM_OPT_MODE                  2\n")
    fsm_file.write("METHOD                        %s\n" %method)
    fsm_file.write("BASIS                         %s\n" %basis)
    fsm_file.write("unrestricted                  true\n")
    fsm_file.write("MEM_STATIC                    256\n")
    fsm_file.write("MAX_SCF_CYCLES                500\n")
    fsm_file.write("MEM_TOTAL                     4096\n")
    fsm_file.write("$end")

parser = argparse.ArgumentParser()
parser.add_argument('-f1', help='.xyz file or directory (reactant geometry)', type=str)
parser.add_argument('-f2', help='.xyz file or directory (product geometry)', type=str)
#parser.add_argument('-j', help='Job Type', type=str,default='fsm') # currently the only option is fsm
parser.add_argument('-c', help='Charge', type=int,default=0) # unlike in qcIn.py, -c and -m are ints here
parser.add_argument('-m', help='Multiplicity', type=int,default=1)
parser.add_argument('-basis', help='Basis Set', type=str,default='6-31+g*')
parser.add_argument('-method', help='EST Method', type=str,default='tpssh')
parser.add_argument('-node', help='number of nodes for FSM', type=int,default=20)
args = parser.parse_args()

f1 = args.f1.rstrip('/')
f2 = args.f2.rstrip('/')

if os.path.isfile(f1) and os.path.isfile(f2):
    if f1.split('.')[-1] == "xyz" and f2.split('.')[-1] == "xyz":
        write_fsm_file(f1, f2, args.c, args.m, args.basis, args.method, args.node)
    else:
        print("Need to input two .xyz files")
elif os.path.isdir(f1) and os.path.isdir(f2):
    for product_file in os.listdir(f2):
        if product_file.split('.')[-1] != "xyz":
            continue # also if this happens you might be about to overwrite something
        reactant_file = find_reactant_file(product_file, f1)
        if reactant_file == None:
            print("Couldn't find reactant file for " + os.path.basename(product_file))
            continue
        write_fsm_file(f1 + '/' + reactant_file, f2 + '/' + product_file, args.c, args.m, args.basis, args.method, args.node)
else: # including if f1 or f2 doesn't exist
    print("Need to input two files or two directories")
