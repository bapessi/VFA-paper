import pandas as pd 
import numpy as np
import scipy.optimize as SO
import scipy.interpolate as SI
from scipy.integrate import ode

def read_data(excel_name):

    df = pd.read_excel(excel_name)
    df = df.dropna(axis=0,how = 'all')
    df['paper'] = df['paper'].fillna(method='ffill')
    df['n_exp'] = df['n_exp'].fillna(method='ffill')
    return df
    
def get_init_params(file_params):
    df_params = pd.read_excel(file_params,index_col='parameter')
    df_params = df_params.drop('unit',axis=1)
    df_params = df_params.to_dict('dict')
    df_params = df_params['value']
    return df_params

def treat_data(df):
    papers = list(set(df['paper']))
    data_dic = {}
    v_C = {'X': 1, 'GLU': 0.4, 'GLY': 0.39, 'ACE':0.4 ,'BUT': 0.54, 'GAP': 1, 'SUC': 1}
    
    for paper in papers:
        exp_dic = {}
        for exp in set(df.loc[df['paper']==paper]['n_exp']):

            t_points = []
            y0 = np.zeros(len(variables))
            # y0[-1] = 8.178125/(1000*32) This is necessary in the case of Oxygen
            df2 = df.loc[df['paper']==paper]
            df2 = df2.loc[df2['n_exp']==exp]
            df2 = df2.drop(['unit'],axis = 1)
            var_dic = {}
            for v in df2['var']:
                if v =='I':

                    df3 = df2.loc[df2['var']==v]
                    df3 = df3.drop([df3.columns[0],df3.columns[1],df3.columns[2]],axis=1)
                    df3 = np.array(df3.dropna(axis=1,how='all'),dtype=np.float64)[0] 
                    df4 = df2.loc[df2['var']=='time_'+v]
                    df4 = df4.drop([df4.columns[0],df4.columns[1],df4.columns[2]],axis=1)
                    df4 = np.array(df4.dropna(axis=1,how='all'),dtype=np.float64)[0]*24

                    df3 = np.concatenate((df3,df3[-1:]),0)
                    df4 = np.concatenate((df4,df4[-1:]),0)
                    
                    I_f = SI.interp1d(df4,df3,kind="previous",fill_value='extrapolate')

                if v in variables:
                    
                    y0[variables.index(v)] = df2.loc[df2['var']==v]['t0']/v_M[v]/v_C[v] # gC/L to mol/L
                    df3 = df2.loc[df2['var']==v]
                    df3 = df3.drop([df3.columns[0],df3.columns[1],df3.columns[2]],axis=1)
                    df3 = np.array(df3.dropna(axis=1,how='all'))[0]/v_M[v]/v_C[v]
                    df4 = df2.loc[df2['var']=='time_'+v]
                    df4 = df4.drop([df4.columns[0],df4.columns[1],df4.columns[2]],axis=1)
                    df4 = np.array(df4.dropna(axis=1,how='all'))[0]*24

                    var_dic[v] = {'C': df3,'t': df4}
                else:
                    t_points += list(np.array(np.array(df2.loc[df2['var']==v])[0][3:],
                                        dtype='float64'))
            t_points = np.unique(t_points)*24 # *24 from day to hour
            t_points = t_points[~np.isnan(t_points)]
            exp_dic[exp] = {'t_points': t_points,'var': var_dic,'y0': y0, 'I_f': I_f}

        data_dic[paper] = exp_dic

    return data_dic

def integrate_stiff(y0,t_points,pr,I_f):

    solver = ode(ODE_sys)
    solver.set_integrator('lsoda')
    solver.set_f_params(pr,I_f)
    solver.set_initial_value(y0,t_points[0])

    tf = t_points[-1]
    tvals = []
    y = []
    dt = 1/60
    i = 0

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if solver.t >= t_points[i]:
            tvals.append(t_points[i])
            i += 1
            y.append(res)
            
    y = np.array(y).transpose()
    return y,tvals

def integrate_stiff_taylor(y0,t_points,pr,I_f):

    solver = ode(ODE_sys)
    solver.set_integrator('lsoda')
    solver.set_f_params(pr,I_f)
    solver.set_initial_value(y0,t_points[0])

    tf = t_points[-1]
    tvals = []
    y = []
    dt = 1/60
    i = 0

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if solver.t >= t_points[i]:
            i += 1
            y.append(res)
            
    y = np.array(y).transpose()
    return y

def ODE_sys(t,y,pr,I_f):
    # I0= 1
    I0 = I_f(t)

    BUT, ACE, SUC, B, GAP, GLU, GLY = np.zeros(7)
    if y[0] >0 : BUT = y[0]
    if y[1] >0 : ACE = y[1]
    if y[2] >0 : SUC = y[2]
    if y[3] >0 : B = y[3]
    if y[4] >0 : GAP = y[4]
    if y[5] >0 : GLU = y[5]
    if y[6] >0 : GLY = y[6]

    alpha_1 = pr['k_1']*ACE/(pr['KS_1']+ACE)
    alpha_2 = pr['k_2']*BUT/(BUT+pr['k_2']/pr['beta_2']*(BUT/pr['S_opt']-1)**2)*pr['k_d']/(ACE+pr['k_d'])
    alpha_3 = pr['k_3']*I0*(1-np.exp(-pr['beta_3']*B))/(pr['beta_3']*B)

    # alpha_4 = pr['k_4']*GLY/(pr['KS_4']+GLY)
    alpha_4 = pr['k_4']*GLY/(GLY+pr['k_4']/pr['al_4']*(GLY/pr['KS_4']-1)**2)
    alpha_5 = pr['k_5']*GAP
    alpha_6 = pr['k_6']*SUC
    # alpha_7 = pr['k_7']*GLU/(GLU+pr['KS_4'])
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)

    dBUT = -B*alpha_2
    dACE = -B*2*alpha_1
    dSUC = B*(alpha_1 + alpha_2 - 4.08788*alpha_6)

    dGLY = -B*alpha_4
    dGLU = -B*2.07298*alpha_7
    dGAP = B*(alpha_4 -4.08788*alpha_5 + alpha_3)
    dB = B*(alpha_5 + alpha_7 + alpha_6)

    return [dBUT,dACE,dSUC,dB,dGAP,dGLU,dGLY]

variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
