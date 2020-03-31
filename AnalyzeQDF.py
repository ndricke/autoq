import pandas as pd
import sys
import subprocess


def get_restarts(qdf):
    qdf_fail = qdf[(qdf.Esolv.isnull()) | (qdf.GeometryConverged == False)]
    qdf_success = qdf[~((qdf.Esolv.isnull()) | (qdf.GeometryConverged == False))]

    # Find the lowest energy successful job matching the species and charge
    for index, row in qdf_fail.iterrows():
        qdf_spec = qdf_success[qdf_success["Species"] == row["Species"]]
        if qdf_spec.empty:
            print(row["Species"], " has no corresponding successful calculations for taking geometry")
        else:
            qdf_match = qdf_spec[qdf_spec["Charge"] == row["Charge"]]
            if qdf_match.empty:
                qdf_match = qdf_spec
            match_index = qdf_match["Esolv"].idxmin(axis=0)
            rr_job = qdf_match.loc[match_index]
            rr_initial = rr_job["Filename"]
            rr_chmult = row["Filename"].split('.')[0].split('_')[-1]
            rr_final = "%s%s_rr_%s.xyz" % (working_dir+'/', row["Species"], rr_chmult)
            print(working_dir+'/'+rr_initial, rr_final)
            subprocess.run(["babel", "-iqcout", working_dir+'/'+rr_initial, "-oxyz", rr_final]) #, shell=True)


        


infile = sys.argv[1]
working_dir = sys.argv[2]

df = pd.read_json(infile)
get_restarts(df)
