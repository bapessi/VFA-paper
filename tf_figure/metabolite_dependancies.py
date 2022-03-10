import numpy as np 
import pandas as pd

filename = 'bio_noAA'
filename = 'stoich_noexchange'
file_input = filename + '_corrected.xlsx'
output = 'dependencies_' + file_input 

metabolites_in = ['ACE','BUT','GLYC_ext','GLC','ADP','ATP','H','H2O','Pi','NADPH','NADP','cNADP','cNADPH',
    'CO2','O2','SO4','NH4','Mg2','NAD','NADH','ACP','Light','THF','CoA','FAD']
 
df = pd.read_excel(file_input,index_col=0)
df = df.drop('reversible',axis=0)
metabolites = df.index.values
dict = pd.DataFrame(columns=['pos'],index=metabolites)
for m in metabolites:
    if m in metabolites_in:
        dict.loc[m] = 0
    else:
        dict.loc[m] = -100


to_update = True
while to_update: 
    to_update = False 
    for r in df.columns:
        reacts = df.loc[:,r].loc[df[r]<0]
        products = df.loc[:,r].loc[df[r]>0]
        pos_prod_list = []
        pos_reac_list = []
        for p in products.index: 
            pos_prod_list.append(dict.loc[p].values[0])
        for s in reacts.index:
            pos_reac_list.append(dict.loc[s].values[0])


        if np.all([p>=0 for p in pos_prod_list]):
            for s in reacts.index:
                if dict.loc[s].values[0]<0:
                    dict.loc[s] = max(pos_prod_list)+1
                    to_update = True

        if np.all([s>=0 for s in pos_reac_list]):
            for p in products.index:
                if dict.loc[p].values[0]<0:
                    dict.loc[p] = max(pos_reac_list)+1
                    to_update = True

print(dict.loc[dict['pos']>0])

# print(dict.loc[dict['pos']<0])
