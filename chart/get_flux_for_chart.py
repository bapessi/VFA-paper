import pandas as pd
import numpy as np

fileout = 'data/chart_reactions_positions_y1.xlsx'
fileinput = 'data/FBA_main_biomass_y1.xlsx'
df_in = pd.read_excel(fileinput,index_col=0)
df_out = pd.read_excel(fileout,index_col=0)
flux_list = []
for reaction in df_out['Reaction'].values:
    flux = df_in['Width'].loc[reaction]
    flux_list.append(flux)
    
# flux_series = pd.Series(flux_list,name='flux')
df_out.insert(2,"width",np.array(flux_list))
df_out.to_excel(fileout)

