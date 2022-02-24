import efm_chart as ec
import get_stoichiometry as gs 
import pandas as pd
import numpy as np
# gs.save_stoich_excel(filename)
# filename = 'chloroplast.xlsx'
# efms = ec.get_efms(filename)

def get_feasible(df):
    list_not_feasible =  []
    for m in df.index.values:
        if not np.any(df.loc[m].values[df.loc[m].values>0]):
            list_not_feasible.append(m)
        if not np.any(df.loc[m].values[df.loc[m].values<0]):
            list_not_feasible.append(m)
    list_not_feasible.remove('reversible')
    metabolites_ignore = ['ADP','ATP','H','H2O','Pi','NADPH','NADP','cNADP','cNADPH']
    bool = [m in df.index.values for m in metabolites_ignore]
    print(metabolites_ignore)
    # list_not_feasible = list_not_feasible + metabolites_ignore 
    return list_not_feasible 
def del_non_feasible(df,list_not_feasible):
    for m in list_not_feasible:
        df = df.drop(m)
    for r in df.columns.values:
        if not np.any(df.loc[:,r].values[df.loc[:,r].values>0]):
            print('reaction problem')
            df = df.drop(r,axis=1)

    return df 

if __name__ == '__main__':
    filename = 'photo_brut_corrected.xlsx'
    df = pd.read_excel(filename,index_col=0)
    list_not_feasible = get_feasible(df)
    df = del_non_feasible(df,list_not_feasible)
    df.to_excel('photo_feasible.xlsx')
    filename='photo_feasible.xlsx'
    # efms= ec.get_efms(filename)    
    # print(efms)
   
