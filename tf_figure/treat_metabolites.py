# This file gets the excel from word format and remove metabolites and reactions to make feasibles efms calculations
import efm_chart as ec
import get_stoichiometry as gs 
import pandas as pd
import numpy as np

def get_feasible(df,metabolites_ignore):
    list_not_feasible =  []
    for m in df.index.values:
        if not np.any(df.loc[m].values[df.loc[m].values>0]):
            list_not_feasible.append(m)
            print(m,' is not feasible')
        if not np.any(df.loc[m].values[df.loc[m].values<0]):
            list_not_feasible.append(m)
            print(m,' is not feasible')
    list_not_feasible.remove('reversible')
    bool = [m in df.index.values for m in metabolites_ignore]
    met_new = []
    for m,b in zip(metabolites_ignore,bool):
        if b and not m in list_not_feasible:
            met_new.append(m)  
    list_not_feasible = list_not_feasible + met_new 
    return list_not_feasible 

def del_non_feasible(df,list_not_feasible):

    for m in list_not_feasible:
        df = df.drop(m)
    for r in df.columns.values:
        if not np.any(df.loc[:,r].values[df.loc[:,r].values!=0]):
            print('reaction problem')
            df = df.drop(r,axis=1)

    return df 

if __name__ == '__main__':

    filename = 'stoich_biomass_halffull' # Excel directed from word file
    filename = 'TCC'
    filename = 'PPP'
    filename = 'nucleic_acids'
    filename = 'chlorophyll'
    filename = 'THF'
    filename = 'bio_noAA'
    filename = 'AA'
    filename = 'stoich_noexchange'
    filename = 'main_biomass'
    # filename = 'lipid_synthesis' 
    # filename = 'AA100'

    # metabolites_ignore = ['ADP','ATP','H','H2O','Pi','NADPH','NADP','cNADP','cNADPH',
    # 'CO2','O2','SO4','NH4','Mg2','NAD','NADH','ACP']
    # metabolites_ignore = ['BUT','GLYC_ext','GLC','ADP','ATP','H','H2O','Pi','NADPH','NADP','cNADP','cNADPH',
    # 'CO2','O2','SO4','NH4','Mg2','NAD','NADH','ACP','Light','THF','CoA','FAD']
    metabolites_ignore = ['BUT','GLYC_ext','GLC','ADP','ATP','H','H2O','Pi','NADPH','NADP','cNADP','cNADPH',
    'CO2','O2','SO4','NH4','Mg2','NAD','NADH','ACP','Light','CoA','FAD']
 
    file_input = filename + '_corrected.xlsx' # excel format for FBA and EFMs 
    file_output = filename +'_feasible.xlsx' # final excel with feasible reactions and metabolites for FBA and EFMs 
    gs.save_sto_corrected(filename) # treat excel brut from word
    df = pd.read_excel(file_input,index_col=0)
    list_not_feasible = get_feasible(df,metabolites_ignore)
    df = del_non_feasible(df,list_not_feasible)
    df.to_excel(file_output)

    efms= ec.get_efms(file_output)    
    # print(np.column_stack((efms,df.columns.values)))
    print(efms)
    print("number of emfs: ", len(efms.T))
