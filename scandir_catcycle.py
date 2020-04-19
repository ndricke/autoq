import sys
import pandas as pd
import scandir_analysis


df = pd.read_json(sys.argv[1])
df_bound_bare = scandir_analysis.calc_binding_energy_autoq(df)
df_bound_bare.to_json(sys.argv[2])


