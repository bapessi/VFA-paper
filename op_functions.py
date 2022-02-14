import numpy as np
from scipy.integrate import ode
import scipy.optimize as SO
I0 = 500/136
def ODE_sys(t,y,pr):
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


def min_stiff(y0,t_points,pr):


    variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
    solver = ode(ODE_sys)
    solver.set_integrator('lsoda')
    solver.set_f_params(pr)
    solver.set_initial_value(y0,t_points[0])
    tf = t_points[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <=0.1027/92 :

            break

    
    return solver.t

def min_time(g0,v_M,pr):
    pr = pr
    GLU0,GLY0 = g0 
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]

    t_points = [0,50*24]
    tf = min_stiff(y0,t_points,pr)
    #print('GLU0 ', GLU0, ' GLY0 ', GLY0)
    #print('time ',tf/24)     
    return tf/24
