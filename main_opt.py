import sys
sys.path.append('../')
import notebook_functions as nf
import scipy.interpolate as SI
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import scipy.optimize as SO
import pickle
import op_functions as op
import op_fedbatch2 as opf
import op_ta as ota

def get_diff_conditions_batch():
    bounds_list = [[(0,500),(0,500)],[(0,0),(0,500)],[(0,500),(0,0)]] 
    for bounds in bounds_list:
        res = SO.differential_evolution(op.min_time,bounds,workers=-1,args=[v_M,pr])
        print(res)

def get_diff_conditions_fedbatch(substrate):
    print('fed batch')
    if substrate =='GLU':
        fr = 1
    if substrate =='GLY':
        fr = 0
    
    t_points = [0,10*24]
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']
    Ta = np.array([0.5])
    D0 = np.array([0.21459823,1.57187414e-02, 0.02659192, 0.21459823])
    bounds = tuple([(0.1,3)]+len(D0)*[(0,1)])
    params0 = np.concatenate((Ta,D0))
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,0,GLY0/v_M['GLY'],GLU0/v_M['GLU']]

    res = SO.differential_evolution(opf.min_time,bounds=bounds,workers=-1,x0=params0,args=[y0,t_points,pr,fr,S0])
    row =np.array([res.x[0],fr,res.fun,*res.x[1:],S0]).reshape(1,len(res.x)+3)
    return row

def get_diff_fraction_fedbatch():
    print('fed batch')
    t_points = [0,10*24]
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']
    Ta = np.array([0.5])
    fr = np.array([0.5])
    D0 = np.array([0.21459823,1.57187414e-02, 0.02659192, 0.21459823])
    
    bounds = tuple([(0.1,3)]+[(0,1)]+len(D0)*[(0,1)])
    params0 = np.concatenate((Ta,fr,D0))
    print(params0)
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,0,GLY0/v_M['GLY'],GLU0/v_M['GLU']]

    res = SO.differential_evolution(opf.min_time_fraction,bounds=bounds,workers=-1,x0=params0,args=[y0,t_points,pr,S0])
    row =np.array([res.x[1],res.x[0],res.fun,*res.x[2:],S0]).reshape(1,len(res.x)+2)
    return row
def save_data(row):
    df = pd.read_excel(file_results)
    columns = ['Ta','fr','tf','D0','D1','D2','D3','S0'] 
    df_new_results = pd.DataFrame(row,columns=columns)
    df_new_results =pd.concat([df,df_new_results],ignore_index=True)
    df_new_results.to_excel(file_results,index=False)


if __name__ == '__main__':
    file_params = '../Parameters_inhibition.xlsx'
    file_results = 'results_op_fraction.xlsx'
    pr = nf.get_init_params(file_params)
    GLU0 = 0
    GLY0 = 0
    v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    
   # get_diff_conditions_batch()
    S0 = 2000
    row = ota.get_diff_conditions_fedbatch('GLY',pr,S0)
    ota.save_data(row)
    row = get_diff_conditions_fedbatch(substrate='GLU')
    save_data(row)
    row = get_diff_conditions_fedbatch(substrate='GLY')
    save_data(row) 
    row = get_diff_fraction_fedbatch()
    save_data(row)
