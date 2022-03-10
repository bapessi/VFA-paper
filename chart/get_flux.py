import notebook_functions as nf
import numpy as np
import pandas as pd
def get_flux_zero(y):
    file_params = 'data/' + 'Parameters_inhibition.xlsx'
    I0 = 500/136 
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
    alpha_list =[alpha_1,alpha_2,alpha_3,alpha_4,alpha_5,alpha_6,alpha_7]
    return dt,alpha_list

v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
GLY0 = 0.071*v_M['GLY']
GLU0 = 0.051*v_M['GLU']
y = [3.5/v_M['BUT'],1.7/v_M["ACE"],0,0.1/v_M['X'],0,GLU0/v_M['GLU'],GLY0/v_M['GLY']]  
y= [3.80506107e-02, 0,  2.36786062e-08,  9.59648585e-02,
  5.07225032e-05,  5.42885076e-01,  1.85227299e-04]

dt,alpha_list = get_flux_zero(y)
print(dt)
df = pd.DataFrame(np.array([dt,alpha_list]).T,columns=['dt','alpha_list'])
print(df)
df.to_excel('data/flux_values.xlsx')
