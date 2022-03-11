import pandas as pd 
import numpy as np
import notebook_functions as nf

def get_cat(BUT,ACE,GLU,GLY,I):
    if I >1e-4:
        sub = 'MIX'
    else:
        sub = 'HETE'
    if BUT==1 and ACE==1:
        cat='ACE-BUT'
    elif GLY==1 and GLU==1:
        cat='GLY-GLU'
    elif BUT==1:
        cat='BUT'
    elif ACE==1:
        cat='ACE'
    elif GLU==1:
        cat='GLU'
    elif GLY==1:
        cat='GLY'
    else:
        cat='AUTO'
        sub =''
    return cat,sub

def iterate_papers(df):
    papers = list(df.keys())
    cat_experiments = pd.DataFrame(columns=categories)
    exp = list(df[papers[0]].keys())
    n_exp = 0

    for paper in papers:
        list_hete = np.zeros((len(categories)))
        list_mix = np.zeros((len(categories)))
        for exp in df[paper].keys():
            BUT,ACE,GLU,GLY = np.zeros(4)
            y0 = df[paper][exp]['y0']
            I0 = df[paper][exp]['I_f'](0)

            if y0[variables.index('BUT')] > 1e-4:
                BUT = 1
            if y0[variables.index('ACE')] > 1e-4:
                ACE = 1
            if y0[variables.index('GLY')] > 1e-4:
                GLY = 1
            if y0[variables.index('GLU')] > 1e-4:
                GLU = 1

            cat,sub = get_cat(BUT,ACE,GLU,GLY,I0)
            if sub =='HETE':
                list_hete[categories.index(cat)] = list_hete[categories.index(cat)] + 1 
            else:
                list_mix[categories.index(cat)] = list_mix[categories.index(cat)] + 1 

        df2 = pd.DataFrame([list_hete,list_mix],index= [paper+' HETE',paper+' MIX'],columns=categories)
        cat_experiments = pd.concat([cat_experiments,df2],axis=0)

    return cat_experiments 

if __name__ == '__main__':
    filename = 'dataSimulateAll.xls'
    df = nf.read_data(filename)
    df = nf.treat_data(df)

    categories = ('GLU','GLY','GLY-GLU','ACE','BUT','ACE-BUT','AUTO')
    sub_category = ('MIX','HETE')
    variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
    df = iterate_papers(df)
    # df = df.loc[df['AUTO']>0]
    df.to_excel('experiments_categories.xlsx')
    print(df) 
