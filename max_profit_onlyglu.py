import notebook_functions as nf
import fedbatch_Vcte_glugly as fd
import numpy as np
from scipy.integrate import ode
import scipy.optimize as so
from scipy import integrate

def ODE(t,y,Ta):

    S0 = 1000
    Sin_glu = S0/v_M['GLU']
    Sin_gly = 0 
    Sin = Sin_gly+Sin_glu
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
    # alpha_3 = pr['k_3']*I0*(1-np.exp(-pr['beta_3']*B))/(pr['beta_3']*B)
    alpha_3B = pr['k_3']*I0*(1-np.exp(-pr['beta_3']*B))/(pr['beta_3'])
    alpha_4 = pr['k_4']*GLY/(GLY+pr['k_4']/pr['al_4']*(GLY/pr['KS_4']-1)**2)
    alpha_5 = pr['k_5']*GAP
    alpha_6 = pr['k_6']*SUC
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)

    if t < Ta:
        D = B*2.07298*alpha_7/(Sin_glu)
    else:
        D = 0 

    dACE = -B*2*alpha_1 
    dBUT = -B*alpha_2 
    dSUC = B*(alpha_1 + alpha_2 - 4.08788*alpha_6) 
    dGLY = -B*alpha_4 + D*Sin_gly 
    dGLU = -B*2.07298*alpha_7 + D*Sin_glu 
    # dGAP = B*(alpha_4 -4.08788*alpha_5 + alpha_3) 
    dGAP = B*(alpha_4 -4.08788*alpha_5) + alpha_3B
    dB = B*(alpha_5 + alpha_7 + alpha_6) 

    dt = [dBUT,dACE,dSUC,dB,dGAP,dGLU,dGLY]
    return dt


def integrate_ODE(p_control):
    GLY0 = 0
    GLU0 = 0.051*v_M['GLU']   
    Ta= p_control
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    t_eval = np.arange(0,24*6,0.1)
    solver = ode(ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(Ta)
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
    return solver.t/24,B,GLU,t_list 


def get_total_GLU_fedbatch(t_eval,B,GLU):
    S0 = 1000
    Sin_glu = S0/v_M['GLU']
    alpha_7 = pr['k_7']*GLU/(GLU+pr['k_7']/pr['al_7']*(GLU/pr['KS_7']-1)**2)
    D = B*2.07298*alpha_7/(Sin_glu)
    totalGLU =  integrate.simpson(D*GLU,t_eval)*v_M['GLU']
    return totalGLU #g/L


def get_profit_fedbatch(p_control):
    Ta = p_control
    priceGlu = 410
    GLY0 = 0
    GLU0 = 0.051*v_M['GLU']
    priceBio = 1000
    tf,B,GLU,t_eval = integrate_ODE(p_control) 
    GLYtotal = 0
    GLUtotal = get_total_GLU_fedbatch(t_eval,B,GLU) 
    Xf = B[-1]*v_M['X'] #g/L 
    return -(priceBio*Xf- GLUtotal*priceGlu)/(tf*1e+6)

file_params = 'Parameters_inhibition.xlsx'
pr = nf.get_init_params(file_params)
v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
bounds = [(0,24*4)]

if __name__ == '__main__':
    res = so.differential_evolution(get_profit_fedbatch,bounds=bounds)
    print(res)
    print(res.x)
