import sys
import pandas as pd
import numpy as np
import re

"""
A script to be run on the output file from the python script scandir.py
This script assumes it is being applied to ORR catalysts, and makes some assumptions about how the naming scheme
describes the contents of the output files. It would be good to adapt graph-based analysis in the future, but it works
as a quick hack for now.

"""

mol_dict = {"O2": -150.32854, "O2H": 0, "O": 0, "OH":-75.924145, "CO":-113.3163168, "CN":-92.96919536}
mol_list = ["O2", "O2br", "O2r", "O2H", "O", "OH", "CO", "CN", "None"]
mol_energies = [-150.32854, -150.32854, -150.32854, 0, 0, -75.924145, -113.3163168, -92.96919536, 0]
df_molecules = pd.DataFrame({"Bound":mol_list, "Esolv_bound":mol_energies})


def map_size(species):
    if "X2" in species:
        return "X2"
    elif "X" in species:
        return "X"
    else:
        return "Reg"


def map_autoq_catalysts(filename):
    """
    Take a filename, and separate out the bound species
    """
    s_spec = pd.Series(filename)
    catalysts = ("C4Fe", "N2C2Fe", "N2rC2Fe", "nanFe", "nanFe-functionalized", "porphyrinFe-functionalized", "tetrid", "tetry", "mepyr")
    bound = ("O2", "O2r",  "O2br", "O2H", "OH", "CN", "CO", "O")
    s_spec["Catalyst"] = None
    s_spec["Bound"] = None

    # need to fix problem that species does not always contain the func## (accidentally stored in jobtype by some schemes)
    # it maybe be good to just rename those jobs to better follow some aspects of the naming scheme
    # but we have the examples of:
    # tetrid_func9_optsp_a0m2.out --> tetrid
    # tetridOH_func1_optsp_a0m1.out --> tetridOH
    # mepyr-func36O2H_optsp-14_optsp_a0m1.out --> mepyr-func36O2H
    # it might be better just to build something that can fall back on filename with re under certain conditions?
    for catalyst in catalysts:
        if catalyst in filename:
            #print("found %s in %s" % (catalyst, filename))
            s_spec["Catalyst"] = catalyst

            number_match = re.search("func(\d+)([^-|_]*)-?(\d+)?", filename)
            if number_match:
                s_spec["Funcnum"] = int(number_match.group(1))
                s_spec["Bound"] = number_match.group(2)
                s_spec["Bound_site"] = number_match.group(3)

            if s_spec["Bound"] == None or s_spec["Bound"] == '': # this will happen in the case of tetridOH_func5_optsp_a0m1.out
                trunc_species = filename.split('_')[0].replace(catalyst, '').replace("X2", '').replace("X", '')
                if trunc_species in bound:
                    s_spec["Bound"] = trunc_species
                else:
                    print(trunc_species, catalyst)
                
    return s_spec


def map_catalysts(species):
    """
    Take a filename, and separate out the bound species
    """
    s_spec = pd.Series(species)
    catalysts = ("C4Fe", "N2C2Fe", "N2rC2Fe", "nanFe", "nanFe-functionalized", "porphyrinFe-functionalized", "tetrid", "tetry", "mepyr")
    bound = ("O2", "O2r",  "O2br", "O2H", "OH", "CN", "CO", "O")
    s_spec["Catalyst"] = None
    s_spec["Bound"] = None

    for catalyst in catalysts:
        if catalyst in species:
            s_spec["Catalyst"] = catalyst
            trunc_species = species.replace(catalyst, '').replace("X2", '').replace("X", '')
                    
            number_match = re.search(catalyst+"(\d+)", species)
            if number_match:
                catalyst += number_match.group(1)
                s_spec["Funcnum"] = int(number_match.group(1))

            for react in bound:
                if react == trunc_species:
                    s_spec["Bound"] = react

    return s_spec


def parse_bound_by_name(df):
    """For a list of species tags, figure out what was bound based on the name"""
    return df["Species"].apply(map_catalysts)


def parse_size_by_name(df):
    """If X2, set size X2. If only X, set size X. If neither X2 nor X, set regular"""
    return df["Species"].apply(map_size)

def get_min_row(df):
    min_energy_idx = df["Esolv"].idxmin()
    print(df.shape, min_energy_idx)
    return df.iloc[min_energy_idx]

def find_min_bound(df, column_groups):
    """For each charge state, size, catalyst, and bound species, find the min energy state"""
    df_min_idx = df.groupby(column_groups, sort=False)["Esolv"].idxmin()
    df_min = df.loc[df_min_idx]
    return df_min


def calc_binding_energy_autoq(df):
    """Calculate binding energy from bare to bound for each bound species"""
    # Split by bare and bound catalysts
    df["data_dir"] = df # TODO strip Bound off of data_dir
    df_bare = df[df["Bound"] == ""]
    df_bound = df[df["Bound"] != ""]


    df_bound_bare = df_bound.merge(df_bare[["Catalyst", "Funcnum", "data_dir", "Esolv"]], how="outer", 
                    on=["Catalyst", "Funcnum", "data_dir"], suffixes=("", "_bare"))

    df_bound_bare["E_binding"] = df_bound_bare["Esolv"] - df_bound_bare["Esolv_bare"]
    return df_bound_bare


def calc_binding_energies(df, df_bound_species):
    """Calculate binding energy from bare to bound for each bound species"""
    bound_charge_map = {"O2":0, "O2br":0, "O2r":0, "O2H":0, "O":0, "OH":-1, "CO":0, "CN":-1, "None":0}

    # Split by bare and bound catalysts
    df_bare = df[df["Bound"] == "None"]
    df_bound = df[df["Bound"] != "None"]

    df_bare.rename(columns={"Charge":"Charge_Bare"}, inplace=True)

    df_bound["Charge_Mod"] = df_bound["Bound"].map(bound_charge_map)
    df_bound["Charge_Bare"] = df_bound["Charge"] - df_bound["Charge_Mod"]
    df_bound_bare = df_bound.merge(df_bare[["Charge_Bare", "Size", "Catalyst", "Esolv"]], how="outer", 
                    on=["Charge_Bare", "Size", "Catalyst"], suffixes=("", "_bare"))

    print(df_bound_bare.shape)
    df_bound_bare = df_bound_bare.merge(df_bound_species, how="left", on="Bound")
    print(df_bound_bare.shape)
    return df_bound_bare


def parse_HER_catalysts(df):
    df["Size"] = parse_size_by_name(df)
    df_spec_split = parse_bound_by_name(df)
    df_aug = df.merge(df_spec_split, left_on="Species", right_on=0)
    df_aug_min = find_min_bound(df_aug, column_groups=["Charge", "Size", "Catalyst", "Bound"])
    df_fin = calc_binding_energies(df_aug_min, df_molecules)
    df_fin["Binding_Energy"] = (df_fin["Esolv"] - df_fin["Esolv_bare"] - df_fin["Esolv_bound"])*27.211
    df_fin = df_fin.dropna()
    df_fin.to_csv(outfile + "_bindingE.csv")
    df_aug.to_csv(outfile + "_annotated.csv")

def parse_autoq_catalysts(df):
    df = df[~df["Species"].str.contains("_fq_")]
    df_spec_split = df["Filename"].apply(map_autoq_catalysts)
    #print(type(df_spec_split))
    #print(df_spec_split)
    df_aug = df.merge(df_spec_split, left_on="Filename", right_on=0)
    #df_aug_min = find_min_bound(df_aug, column_groups=["Charge", "Catalyst", "Bound"]) # there is no min bound for autoq
    return df_aug

if __name__ == "__main__":

    read_funcs = {"json": pd.read_json, "csv": pd.read_csv}

    infile = sys.argv[1]
    outfile = sys.argv[2]


    read_func = read_funcs[infile.split('.')[-1]]
    df = read_func(infile)
    #parse_HER_catalysts(df)
    parse_autoq_catalysts(df)














