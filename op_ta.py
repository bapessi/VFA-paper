from scipy import interpolate
import pandas as pd
import numpy as np
import scipy.optimize as SO
from scipy.integrate import ode
v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']

def get_diff_conditions_fedbatch(substrate,pr,S0):
    print('fed batch')
    if substrate =='GLU':
        fr = 1
    if substrate =='GLY':
        fr = 0
    
    t_points = [0,10*24]
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']
    Ta = np.array([0.5])
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLY0/v_M['GLY'],GLU0/v_M['GLU'],1]
    bounds = tuple([(0.1,3)])

    res = SO.differential_evolution(min_time,bounds=bounds,workers=-1,args=[y0,t_points,pr,fr,S0])
    row =np.array([res.x[0],fr,res.fun,*res.x[1:],S0,' ']).reshape(1,len(res.x)+4)
    return row

def min_time(params,y0,t_points,pr,fr,S0):
    Ta = params[0]
    print('Ta = ',Ta)
    #t_vector = 24*np.array([Ta,Ta+0.25,Ta+0.5,Ta+0.75,Ta+1])
    t_vector = 24*np.array([Ta,Ta+0.5,Ta+1.0,Ta+2.0])
    fr = 0.5
    tf = min_stiff(y0,t_points,Ta,pr,fr,S0)
    print('tf = ',tf/24)    
    return tf

def min_stiff(y0,t_points,Ta,pr,fr,S0):
      
    solver = ode(ODE_sys_Ta_fed_batch)
    solver.set_integrator('lsoda')
    solver.set_f_params(Ta,pr,fr,S0)
    solver.set_initial_value(y0,t_points[0])
    tf = t_points[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <=0.1027/92 :
           
            break
    
    return solver.t



def ODE_sys_Ta_fed_batch(t,y,Ta,pr,fr,S0):
    # DS = D_f(t)
    Sin_glu = S0/v_M['GLU']*fr
    Sin_gly = S0/v_M['GLY']*(1-fr)
    Sin = Sin_gly + Sin_glu 
    I0 = 500/136 
    BUT, ACE, SUC, B, GAP, GLU, GLY = np.zeros(7)
    if y[0] >0 : BUT = y[0]
    if y[1] >0 : ACE = y[1]
    if y[2] >0 : SUC = y[2]
    if y[3] >0 : B = y[3]
    if y[4] >0 : GAP = y[4]
    if y[5] >0 : GLU = y[5]
    if y[6] >0 : GLY = y[6]
    V = y[7]
    
    alpha_1 = pr['k_1']*ACE/(pr['KS_1']+ACE)
    alpha_2 = pr['k_2']*BUT/(BUT+pr['k_2']/pr['beta_2']*(BUT/pr['S_opt']-1)**2)*pr['k_d']/(ACE+pr['k_d'])
    alpha_3 = pr['k_3']*I0*(1-np.exp(-pr['beta_3']*B))/(pr['beta_3']*B)

    alpha_4 = pr['k_4']*GLY/(GLY+pr['k_4']/pr['al_4']*(GLY/pr['KS_4']-1)**2)
    alpha_5 = pr['k_5']*GAP
    alpha_6 = pr['k_6']*SUC
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)

    if t < Ta*24:
        D = B*4.08788*alpha_4/(Sin_gly+GLY)
        #D = B*2.07298*alpha_7/(Sin_glu+GLU)

    else:
        D = 0

    dBUT = -B*alpha_2 -D*BUT
    dACE = -B*2*alpha_1 -D*ACE
    dSUC = B*(alpha_1 + alpha_2 - 4.08788*alpha_6) -D*SUC

    dGLY = -B*alpha_4 -D*GLY +D*Sin_gly
    dGLU = -B*2.07298*alpha_7 + D*Sin_glu - D*GLU

    dGAP = B*(alpha_4 -4.08788*alpha_5 + alpha_3) -D*GAP 
    dB = B*(alpha_5 + alpha_7 + alpha_6) -D*B
    dV = D*V
    dt = [dBUT,dACE,dSUC,dB,dGAP,dGLU,dGLY,dV]
    return dt

def save_data(row):
    file_results_Ta = 'results_Ta.xlsx'
    df = pd.read_excel(file_results_Ta)
    columns = ['Ta','fr','tf','S0','comments']
    df_new_results = pd.DataFrame(row,columns=columns)
    df_new_results =pd.concat([df,df_new_results],ignore_index=True)
    df_new_results.to_excel(file_results_Ta,index=False)


