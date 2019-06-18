import pandas as pd
import os
import sys
def getAttributes(attribute, filename):
    with open(filename, "r") as file:
        chelpg_charge_list = []
        total_free_energy = 0
        line = file.readline()
        while line:
            if attribute=="Mulliken":
                keyphrase = "Mulliken"
                breakphrase = "---"
                deathphrase = "Q-Chem fatal error"
            elif attribute=="ChElPG":
                keyphrase = "ChElPG Net"
                breakphrase = "---"
                deathphrase = "Q-Chem fatal error"
            elif attribute=="Energy":
                keyphrase = "Total Free Energy"
                deathphrase = "Q-Chem fatal error"
                breakphrase = "---"
            if deathphrase in line:
                #print("This file has a fatal error! Moving on...")
                return("Move on")
            if keyphrase in line:
                if attribute=="Energy":
                    splitline = line.split()
                    total_free_energy = splitline[9]
                atom_list = []
                charge_list = []
                spin_list = []
                atom_id_list = []
                counter = 0
                for x in range(3):
                    line = file.readline()
                while True:
                    counter+=1
                    line = file.readline()
                    if breakphrase in line:
                        break
                    splitline = line.split()
                    if attribute=="Mulliken":
                        atom_id_list.append(counter)
                        atom_list.append(splitline[1])
                        charge_list.append(float(splitline[2]))
                        spin_list.append(float(splitline[3]))                                      
                    elif attribute=="ChElPG":
                        chelpg_charge_list.append(float(splitline[2]))
            line = file.readline()
    if attribute=="Mulliken":
        mulliken_list = [atom_id_list, atom_list, charge_list, spin_list]
        # print("Mulliken data collected!")
        # print(mulliken_list)
        return(mulliken_list)
    elif attribute=="ChElPG":
        #print("ChElPG data collected!")         
        return(chelpg_charge_list)
    elif attribute=="Energy":
        #print("Total free energy collected!")
        return(total_free_energy)

def makeDataFrame(input):
    successes = 0
    failures = 0
    df = pd.DataFrame()
    for file in os.listdir(input):
        #print(file)
        fullfilepath = os.path.join(input, file)
        if file.endswith(".out"):
            # print("Processing", file)
            mulliken_part = getAttributes("Mulliken", fullfilepath)
            if mulliken_part == "Move on":
                failures+=1
                continue
            chelpg_part = getAttributes("ChElPG", fullfilepath)
            mulliken_part.append(chelpg_part)
            total_free_energy = getAttributes("Energy", fullfilepath)
            #print("Generating pandas dataframe...")
            length = len(mulliken_part[0])
            catalyst_list = [file]*length
            free_energy_list = [total_free_energy]*length
            all_data = pd.DataFrame(
                {'Atom ID': mulliken_part[0],
                'Atoms': mulliken_part[1],
                'Mulliken Charge': mulliken_part[2],
                'Mulliken Spin': mulliken_part[3],
                'ChElPG Charge': mulliken_part[4],
                'Total Free Energy': free_energy_list,
                'Catalyst': catalyst_list
                })
            df = df.append(all_data)
            successes+=1
    #print(df)
    totalfiles = successes+failures
    print(totalfiles, " files processed")
    print(successes, " successful processes, and ", failures, "failures")
    print("Dataframe created!")
    #print(df)
    return df
if __name__ == "__main__":
    #print(df)
    #print("Running ", sys.argv[0])
    import argparse
    parser = argparse.ArgumentParser("Grab data about bare catalysts and store in csv")
    parser.add_argument('-cat', help='input file/directory (bare catalyst .out)', type=str)
    parser.add_argument('-out', help='output file/directory', type=str)

    args = parser.parse_args()
    inputdirectory = args.cat
    outputdirectory = args.out
    finaldf = makeDataFrame(inputdirectory)
    finalpath = os.path.join(outputdirectory,  'catalyst_data.csv')
    print("Now writing to csv...")
    finaldf.to_csv(finalpath)
