import numpy as np
import pandas as pd
import notebook_functions as nf
from scipy.optimize import linprog

def get_matrix(filename):
    df = pd.read_excel(filename,index_col=0)
    rev = df.loc['reversible',:].values
    df = df.drop('reversible')
    feed_list = determine_feed_list(df)
    Aeq = df.to_numpy(na_value=0)
    beq = np.zeros(Aeq.shape[0]) #  
    c = np.zeros(Aeq.shape[1])
    return df,Aeq,beq,c,rev,feed_list

def get_bounds(rev):
    bounds = []
    # b_max = 2.215 #y0
    b_max = 0.249 #y1
    for r in rev:
        if r == 1:
            bounds.append((-b_max,b_max))
        elif r==0:
            bounds.append((0,b_max))
        else:
            print("error with reversible array")
    return bounds

def correct_feed(feed_list,bounds,c,r_max):
    for reaction,i in zip(df.columns,range(len(df.columns))):
        if reaction in feed_list:
            bounds[i] = (0,2)
            if reaction =='R36': #SUCg to SUC
                # bounds[i] = (0,0.3789/4.08788) # 0.0927  y0
                bounds[i] = (0,0.0023) # y1 
            if reaction == 'R14': # cGAP to GAP
                # bounds[i] = (0,(1.009+0.5259)/4.08788) # 0.3755
                bounds[i] = (0,0.005684) # y1 photo plus gly
                # bounds[i] = (0,0.004) # y1 for photo only
            if reaction == 'R167': #GLU
                # bounds[i] = (0,0.068) #y0
                bounds[i] = (0,0.048) #y1
            if reaction == 'R15': # GLY

                # bounds[i] = (0,1.009) #y0
                bounds[i] = (0,0.0017) #y1
            if reaction == 'R70': #So4
                bounds[i]=(0,10)
            if reaction == 'R20': #BUT
                # bounds[i]=(0,4e-11) # y0
                bounds[i]=(0,0.0023) # y1
            if reaction == 'R27': # ACE
                # bounds[i]=(0,0.38)  # y0
                bounds[i]=(0,1e-11)  # y1
        if reaction == r_max:
            print('max reaction: ',reaction)
            c[i] = -1
    return bounds,c

def determine_feed_list(df):
    feed_list = []
    for r in df.columns:
        if not np.any(df.loc[:,r].values[df.loc[:,r].values < 0]):
            feed_list.append(r) 
    return feed_list

def feed_excel():
    filename = 'data/feed_list.xlsx'
    feed_list = pd.read_excel(filename).to_numpy()
    feed_list=feed_list.flatten()
    return feed_list

if __name__ == '__main__':
    filename = 'main_biomass' + '_feasible' 
    r_max = 'R220'
    # filename = 'glyoxysome' + '_feasible' 
    # r_max = 'R36'
    # filename = 'photo' + '_feasible' 
    # r_max = 'R14'


    fileinput = 'data/' + filename + '.xlsx' 
    fileoutput = 'data/' + 'FBA_' + filename + '.xlsx'
    df,Aeq,beq,c,rev,feed_list = get_matrix(fileinput)
    feed_list = feed_excel()
    bounds = get_bounds(rev)
    bounds,c = correct_feed(feed_list,bounds,c,r_max)
    res = linprog(c,A_eq=Aeq,b_eq=beq,bounds=bounds,method="interior-point",options={'tol':1e-10,'maxiter':100000})
    print("max flux ", max(res.x))
    print("min flux", min(res.x))
    flux =pd.Series(res.x,index=df.columns)
    flux.to_excel(fileoutput)
    print('Flux biomass ',flux['R220'])
    print('Flux 49 ',flux['R49'])
    print('Flux 43,119,121 ',flux['R43'],flux['R119'],flux['R121'])
    print('Flux GLC,GAP, SUC', flux['R167'],flux['R14'],flux['R36'])
