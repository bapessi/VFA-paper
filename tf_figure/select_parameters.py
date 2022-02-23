import numpy as np
import pandas as pd
import notebook_functions as nf
import batch_optimal as bo
import Vglu
import Vgly
import tf_fedbatch as Vglugly

def get_bounds_acebut(pr,params):
    k1p = (0.975*pr['k_1'],1.025*pr['k_1']) # percentage
    KS1p = (1e-7,2.5*5.5e-5) # real value
    k2p = (0.935*pr['k_2'],1.065*pr['k_2']) # percentage
    beta2p = (0.63*pr['beta_2'],1.28*pr['beta_2']) # percentage,
    Soptp = (0.826*pr['S_opt'],1.065*pr['S_opt'])  #percentage
    bounds = [k1p,KS1p,k2p,beta2p,Soptp]
    return bounds

def return_parameters(pr):
    
    glu_samples = pd.read_table('glu_samples.csv',sep=",")
    gly_samples = pd.read_table('gly_samples',sep=",")
    sl = np.random.randint(len(glu_samples))
    gl = np.random.randint(len(gly_samples))
    params_glu = ['k_7','KS_7','al_7']
    params_gly = ['k_4','KS_4','al_4']
    params_glyglu = ['k','KS','al']

    for p,pp in zip(params_glu,params_glyglu):
        pr[p] = glu_samples.iloc[sl][pp]
    for p,pp in zip(params_gly,params_glyglu):
        pr[p] = gly_samples.iloc[gl][pp]

    params_acebut = ['k_1','KS_1','k_2','beta_2','S_opt']
    bounds_acebut = get_bounds_acebut(pr,params_acebut)
    for p,bound in zip(params_acebut,bounds_acebut):
        low,high = bound 
        pr[p]= np.random.uniform(low=low, high=high)

    return pr

def tf_batch():
    # S0 = (112,18) # minimal tf
    S0 = (0,0) # no additional substrate
    tf_list = np.zeros((size))
    pr_old = nf.get_init_params(file_params)
    for i in range(size):
        pr = return_parameters(pr_old.copy()) 
        tf = bo.integrate_ODE(S0,pr)
        tf_list[i] = tf
    return tf_list


def tf_fedbatch_glu():
    size = 100
    Ta_glu = 3.16*24
    tf_list = np.zeros((size))
    pr_old = nf.get_init_params(file_params)
    for i in range(size):
        pr = return_parameters(pr_old.copy()) 
        tf = Vglu.integrate_ODE(Ta_glu,pr)
        tf_list[i] = tf
    return tf_list


def tf_fedbatch_gly():

    size = 100
    Ta_gly = 1.40*24
    tf_list = np.zeros((size))
    pr_old = nf.get_init_params(file_params)
    for i in range(size):
        pr = return_parameters(pr_old.copy()) 
        tf = Vgly.integrate_ODE(Ta_gly,pr)
        tf_list[i] = tf
    return tf_list

def tf_fedbatch_glugly():
    # S0 = (112,18)
    size = 100
    Ta_gly = 1.40*24
    Ta_glu = 1.14*24
    tf_list = np.zeros((size))
    pr_old = nf.get_init_params(file_params)
    for i in range(size):
        pr = return_parameters(pr_old.copy()) 
        tf = Vglugly.integrate_ODE((Ta_glu,Ta_gly),pr)
        tf_list[i] = tf
    return tf_list


if __name__ == '__main__':

    file_params = "Parameters_inhibition.xlsx"
    pr = nf.get_init_params(file_params)
    pr_old = pr.copy()
    v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
    variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']
    # tf_list = tf_batch()
    tf_list = tf_fedbatch_gly()
    print('GLY fedbatch')
    print(np.average(tf_list)/24)
    print(np.std(tf_list)/24)
    # tf = pd.Series(tf_list)
    # tf.to_excel('tf_batch.xlsx')

