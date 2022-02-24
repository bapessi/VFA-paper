import efmtool
import numpy as np
import pandas as pd

def provost_example():

    filename = 'provost_example.xlsx'
    df = pd.read_excel(filename,index_col=0)
    sto = df.to_numpy()
    rev = [0,0,0,0,0,0,0,0]
    reactions = list(np.arange(1,sto.shape[1]+1))
    metabolites = list(np.arange(1,sto.shape[0]+1))
    for m in metabolites:
        metabolites[metabolites.index(m)] = 'M'+str(m) 
    for r in reactions:
        reactions[reactions.index(r)] = 'R'+str(r)
    return sto,rev,reactions,metabolites

def get_sto(filename):
    df = pd.read_excel(filename,index_col=0)
    rev = df.loc['reversible',:].values
    df = df.drop('reversible')
    sto = df.to_numpy()
    metabolites = df.index.values
    reactions = df.columns.values
    return sto,rev,reactions,metabolites

filename = 'PPP_corrected.xlsx'
# filename = 'provost_example.xlsx'
sto,rev,reactions,metabolites = get_sto(filename)

efms = efmtool.calculate_efms(stoichiometry=sto,reversibilities=rev,reaction_names=reactions,metabolite_names=metabolites)
# print(efms)
