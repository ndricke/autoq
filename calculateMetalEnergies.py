import pandas as pd
import os
import sys
def getIntermediateEnergies(filename):
    with open(filename, "r") as file:
        total_free_energy = 0
        line = file.readline()
        keyphrase = "Total Free Energy"
        deathphrase = "Q-Chem fatal error"
        while line:
            if deathphrase in line:
                #print("This file has a fatal error! Moving on...")
                return("Move on")
            if keyphrase in line:
                splitline = line.split()
                total_free_energy = splitline[9]
            line = file.readline()
        return total_free_energy

def makeDataFrame(indir):
    file_list = []
    energy_list = []
    for file in os.listdir(indir):
        full_path = os.path.join(indir, file)
        energy = getIntermediateEnergies(full_path)
        energy_list.append(energy)
        file_list.append(file)
    df = pd.DataFrame(
        {'File': file_list,
        'Energy': energy_list
    
    
    
    })
    print(df)
    return (df)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Grab data about metal catalyst intermediates and store in csv")
    parser.add_argument('-cat', help='input file/directory (catalyst intermediate .out files)', type=str)
    parser.add_argument('-out', help='output directory', type=str)

    args = parser.parse_args()

    outdir = os.path.join(args.out, 'catalyst_energies.csv')
    alldata = makeDataFrame(args.cat)
    alldata.to_csv(outdir)