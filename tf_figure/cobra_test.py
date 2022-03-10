import cobra
import cobra.test
import numpy as np
import pandas as pd
import notebook_functions as nf
from scipy.optimize import linprog

def get_flux_zero():
    file_params = 'Parameters_inhibition.xlsx'
    v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']
    I0 = 500/136 
    y = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]  
    pr = nf.get_init_params(file_params)
    BUT = y[0]
    ACE = y[1]
    SUC = y[2]
    B = y[3]
    GAP = y[4]
    GLU = y[5]
    GLY = y[6]
    
    alpha_1 = pr['k_1']*ACE/(pr['KS_1']+ACE)
    alpha_2 = pr['k_2']*BUT/(BUT+pr['k_2']/pr['beta_2']*(BUT/pr['S_opt']-1)**2)*pr['k_d']/(ACE+pr['k_d'])
    alpha_3 = pr['k_3']*I0*(1-np.exp(-pr['beta_3']*B))/(pr['beta_3']*B)
    alpha_4 = pr['k_4']*GLY/(GLY+pr['k_4']/pr['al_4']*(GLY/pr['KS_4']-1)**2)
    alpha_5 = pr['k_5']*GAP
    alpha_6 = pr['k_6']*SUC
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)

    dACE = -B*2*alpha_1 
    dBUT = -B*alpha_2 
    dSUC = B*(alpha_1 + alpha_2 - 4.08788*alpha_6) 
    dGLY = -B*alpha_4 
    dGLU = -B*2.07298*alpha_7  
    dGAP = B*(alpha_4 -4.08788*alpha_5 + alpha_3) 
    dB = B*(alpha_5 + alpha_7 + alpha_6) 
    dt = [dBUT,dACE,dSUC,dB,dGAP,dGLU,dGLY]

    return dt

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
    for r in rev:
        if r == 1:
            bounds.append((None,None))
        elif r==0:
            bounds.append((0,None))
        else:
            print("error with reversible array")
    return bounds

def correct_feed(feed_list,bounds,c,r_max):
    for reaction,i in zip(df.columns,range(len(df.columns))):
        if reaction in feed_list:
            print('Feed reaction: ',reaction)
            bounds[i] = (0,1e3)
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

if __name__ == '__main__':
    filename = 'glycolysis' + '_feasible' 
    filename = 'TCC' + '_feasible' 
    filename = 'PPP' + '_feasible' 
    filename = 'AA100' + '_feasible' 
    r_max = 'R110' 
    filename = 'nucleic_acids' + '_feasible' 
    r_max = 'R156' 
    filename = 'chlorophyll' + '_feasible' 
    r_max = 'R163' 
    filename = 'THF' + '_feasible' 
    r_max = 'R125' 
    filename = 'bio_noAA' + '_feasible' 
    r_max = 'R166'
    filename = 'lipid_synthesis' + '_feasible' 
    r_max = 'R139'
    filename = 'stoich_noexchange' + '_feasible' 
    r_max = 'R166'
    filename = 'main_biomass' + '_feasible' 
    r_max = 'R166'
    # filename = 'lipid_synthesis' + '_feasible' 
    # filename = 'AA' + '_feasible' 
    # r_max = 'R110' 

    fileinput = filename + '.xlsx' 
    fileoutput = 'FBA_' + filename + '.xlsx'
    df,Aeq,beq,c,rev,feed_list = get_matrix(fileinput)
    bounds = get_bounds(rev)
    bounds,c = correct_feed(feed_list,bounds,c,r_max)
    res = linprog(c,A_eq=Aeq,b_eq=beq,bounds=bounds,method="interior-point",options={'tol':1e-10,'maxiter':10000})
    print(res)
    print("max flux ", max(res.x))
    print("min flux", min(res.x))
    pd.Series(res.x,index=df.columns).to_excel(fileoutput)
    print('feed list', feed_list)
