import Vglu
import notebook_functions as nf
import fedbatch_Vcte_glugly as fd
import numpy as np
from scipy.integrate import ode
import scipy.optimize as so
from scipy import integrate

def get_profit_fedbatch(p_control):
    Ta,fr = p_control
    priceGly = 170
    priceGlu = 410
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']
    priceBio = 1000
    tf,B,GLU,GLY,t_eval = integrate_ODE(p_control) 
    GLYtotal = get_total_GLY_fedbatch(t_eval,B,GLY,fr) 
    GLUtotal = get_total_GLU_fedbatch(t_eval,B,GLU,fr) 
    Xf = B[-1]*v_M['X'] #g/L 
    return -(priceBio*Xf- GLUtotal*priceGlu - GLYtotal*priceGly)/(tf*1e+6)

def integrate_ODE(p_control):
    GLY0 = 0.071*v_M['GLY']
    GLU0 = 0.051*v_M['GLU']   
    Ta,fr = p_control
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    t_eval = np.arange(0,24*6,0.1)
    solver = ode(fd.ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(Ta,fr)
    solver.set_initial_value(y0, t_eval[0])
    tf = t_eval[-1]
    dt = 1/60
    i = 0
    y=[]
    t_list = []
    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if solver.t >= t_eval[i]:
            i = i+1 
            y.append(res)
            t_list.append(solver.t)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <= 0.1027/92:
            break
    y = np.array(y).transpose()
    B= y[variables.index('X')] 
    GLU= y[variables.index('GLU')] 
    GLY= y[variables.index('GLY')] 
    return solver.t/24,B,GLU,GLY,t_list 

def get_total_GLY_fedbatch(t_eval,B,GLY,fr):
    S0 = 1000
    Sin_gly = S0/v_M['GLY']*(1-fr)
    alpha_4 = pr['k_4']*GLY/(GLY+pr['k_4']/pr['al_4']*(GLY/pr['KS_4']-1)**2)
    D = B*alpha_4/Sin_gly
    totalGLY =  integrate.simpson(D*GLY,t_eval)*v_M['GLY']
    return totalGLY

def get_total_GLU_fedbatch(t_eval,B,GLU,fr):
    S0 = 1000
    Sin_glu = S0/v_M['GLU']*fr
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)
    D = B*2.07298*alpha_7/(Sin_glu)
    totalGLU =  integrate.simpson(D*GLU,t_eval)*v_M['GLU']
    return totalGLU #g/L

if __name__ == 'main':
    file_params = 'Parameters_inhibition.xlsx'
    pr = nf.get_init_params(file_params)
    v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
    variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
    bounds = [(0,24*4),(0,1)]
    res = so.differential_evolution(get_profit_fedbatch,bounds=bounds)
    print(res)
    print(res.x)
