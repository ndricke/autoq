import sys
import os
import pandas as pd

import ChemData as CD
import Qdata as QD

indir = sys.argv[1]
outfile = sys.argv[2]

def generate_done(qdatas):
    done_list = []
    not_done_list = []
    for key, item in data.items():
        if item.geo_converged and item.Esolv != None:
            done_list.append(key)
        else:
            not_done_list.append(key)
    return done_list, not_done_list

def generate_qdata_df(qdatas):
    """
    Convert a dictionary of qdatas into a pandas dataframe
    Input
    qdatas (dictionary): filename --> qdata object
    Output
    (pandas dataframe): select data contained in qdatas
    """

    species_list, jobtype_list, charge_list, mult_list, Esolv_list, geo_list, filename_list = [], [], [], [], [], [], []
    for key, qdata in qdatas.items():
        spl_key = key.split('.')[0].split('_')
        species, jobtype = spl_key[0], spl_key[1]
        species_list.append(species)
        jobtype_list.append(jobtype)
        filename_list.append(key)

        # Esolv, geo_converged
        Esolv_list.append(qdata.Esolv)
        geo_list.append(qdata.geo_converged)
        charge_list.append(qdata.charge)
        mult_list.append(qdata.mult)
    

    return pd.DataFrame({"Species":species_list, "JobType":jobtype_list, "Charge":charge_list, "Multiplicity":mult_list, 
                         "Filename":filename_list,
                         "Esolv":Esolv_list, "GeometryConverged":geo_list})


if __name__ == "__main__":
    data = {}
    for filename in os.listdir(indir):
        qcout = indir+'/'+filename
        if qcout.split('.')[-1] == "out":
            qdata = QD.Qdata()
            qdata.readFile(qcout)
            data[filename] = qdata



    df = generate_qdata_df(data)
    df.to_json(outfile+".json")
    df.to_csv(outfile+".csv")

    
