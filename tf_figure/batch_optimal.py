import notebook_functions as nf
import scipy.optimize as so
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import ode

def ODE(t,y,pr):

    I0 = 500/136 

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

def integrate_ODE_Xf(S0):
    GLU0,GLY0 = S0

    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    t_eval = np.linspace(0,24*10,100)
    solver = ode(ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(pr)
    solver.set_initial_value(y0, t_eval[0])
    tf = t_eval[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <= 0.1027/92:
            break

    return solver.t,res[variables.index('X')]*v_M['X']


def integrate_ODE(S0,pr):
    GLU0,GLY0 = S0

    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    t_eval = np.linspace(0,24*30,100)
    solver = ode(ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(pr)
    solver.set_initial_value(y0, t_eval[0])
    tf = t_eval[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <= 0.1027/92:
            break

    return solver.t 

def get_batch_optimal_concentration():
    bounds= [(0,300),(0,300)]
    res = so.differential_evolution(integrate_ODE,bounds=bounds,args=[pr])
    
    return res

def get_profit(S0):
    GLU0,GLY0 = S0
    priceGly = 170 
    priceGlu = 410
    priceBio = 1000 # US$/ton
    tf,Xf = integrate_ODE_Xf(S0)
    profit = (priceBio*Xf- GLU0*priceGlu - GLY0*priceGly)/(tf*1e+6)
    return -profit

def get_batch_optimal_profit():
    bounds= [(0,1000),(0,1000)]
    res = so.differential_evolution(get_profit,bounds=bounds)
    return res


v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']

if __name__ == '__main__':
    file_params = "Parameters_inhibition.xlsx"
    pr = nf.get_init_params(file_params)
    res = get_batch_optimal_concentration()
    print('concentration min tf: ',res.x)
    print("time min tf : ",integrate_ODE(res.x,pr))
    res = get_batch_optimal_profit()
    print("concentration best profit: ",res.x)
    print("time best profit: ",integrate_ODE(res.x,pr))
    print("batch no substrat ", integrate_ODE((0,0),pr))
