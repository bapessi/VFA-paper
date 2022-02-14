from scipy.integrate import ode
import numpy as np
import pandas as pd
import notebook_functions as nf
import matplotlib.pyplot as plt
import scipy.optimize as so
from scipy import interpolate
def ODE(t,y,Ta,D_f):
    S0 = 1000
    Sin_glu = S0/v_M['GLU']
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

    if t < Ta:
        D = B*2.07298*alpha_7/(Sin_glu)
    else:
        D = D_f(t)/S0 

    dACE = -B*2*alpha_1  
    dBUT = -B*alpha_2 
    dSUC = B*(alpha_1 + alpha_2 - 4.08788*alpha_6) 
    dGLY = -B*alpha_4 
    dGLU = -B*2.07298*alpha_7 +D*Sin_glu
    dGAP = B*(alpha_4 -4.08788*alpha_5 + alpha_3) 
    dB = B*(alpha_5 + alpha_7 + alpha_6)
    dt = [dBUT,dACE,dSUC,dB,dGAP,dGLU,dGLY]
    return dt

def integrate_ODE(params):
    Ta = params[0]
    D_array = params[1:] 
    t_vector = np.array([Ta,Ta+5,Ta+8,Ta+10])
    D_f = interpolate.interp1d(t_vector, D_array,kind='linear',fill_value=np.array([0]),bounds_error=False) 
    solver = ode(ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(Ta,D_f)
    solver.set_initial_value(y0, t_eval[0])
    tf = t_eval[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <= 0.1027/92:
            break
    return solver.t 

if __name__ == '__main__':
    file_params = 'Parameters_inhibition.xlsx'
    pr = nf.get_init_params(file_params)
    v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
    variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']

    Ta = np.array([24])
    D0 = np.array([0.21459823,1.57187414e-02, 0.02659192,0.2])
    bounds = tuple([(10,4*24)]+len(D0)*[(0,20)])
    params0 = np.concatenate((Ta,D0))

    GLY0 = 0
    GLU0 = 0.051*v_M['GLU']
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    t_eval = np.linspace(0,10*24,100)
    res = so.differential_evolution(integrate_ODE,bounds=bounds, x0=params0,workers=-1)
    print(res.x)
    Ta = res.x[0]
    tf = integrate_ODE(res.x)
    print('tf = ', tf/24 )
# plt.show()
print('Ta:', Ta/24)
