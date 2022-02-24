import pandas as pd
import numpy as np

filename = 'PPP' 

file_input = filename +'.xlsx'
file_output = filename +'_corrected.xlsx'
df = pd.read_excel(file_input,index_col=0,header=None)

def get_metabolites(df):

    matrix = df.to_numpy()
    metabolites = []
    for line in matrix:
        for cell in line:
            if type(cell)== str:
                if not cell =='+' and not cell=='-->' and not cell =='<-->':
                    metabolites.append(cell)
    metabolites = list(set(metabolites))
    metabolites.sort()
    stoich = pd.DataFrame(index=metabolites,columns=df.index)
    return metabolites,stoich

metabolites,stoich = get_metabolites(df)
reac_rev = np.zeros(len(df.index))
for row in df.index: # df.index are the reactions
    prod = -1
    s = 1
    for col in df.columns: #different points in the initial excel
        if df.at[row,col] in metabolites:
            stoich.at[df.at[row,col],row] = prod*s
        elif df.at[row,col]=='-->' or  df.at[row,col]=='<-->':
            prod = 1
            s = 1
            if df.at[row,col]=='<-->':
                reac_rev[list(df.index).index(row)] = 1

        if df.at[row,col]=='+':
            s = 1
        if type(df.at[row,col])== float or type(df.at[row,col])==int or type(df.at[row,col])==np.float64:
            if not np.isnan(s):
                s = df.at[row,col]

stoich = stoich.fillna(0)
reac_rev = pd.Series(reac_rev,index=df.index,name='reversible')
stoich = stoich.append(reac_rev)
stoich.to_excel(file_output)
print("FINISHED")
