import Vglu
import scipy.optimize as so

def get_profit(x):
    priceGly = 170
    priceGlu = 410
    priceBio = 1000
       
    GLU0, GLY0 = x # g/L
    y0 = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]
    
    tf, Xf = profit_stiff(y0,t_points) # return time (days), X (g/L)
    return -(priceBio*Xf- GLU0*priceGlu - GLY0*priceGly)/(tf*1e+6)

def integrate_ODE(Ta):
    solver = ode(ODE)
    solver.set_integrator('lsoda')
    solver.set_f_params(Ta)
    solver.set_initial_value(y0, t_eval[0])
    tf = t_eval[-1]
    dt = 1/60

    while solver.successful() and solver.t <= tf:
        res = solver.integrate(solver.t+dt)
        if res[variables.index('BUT')] <= 0.069/88 and res[variables.index('GLU')] <= 0.117/180 and res[variables.index('GLY')] <= 0.1027/92:
            break
    return solver.t,res[variables.index('X')]

res = so.differential_evolution(get_profit)
