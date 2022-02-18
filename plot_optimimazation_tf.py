import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import Vgly
import Vglu

glu_samples = pd.read_table('glu_samples_0130.csv',sep=",")
gly_samples = pd.read_table('gly_samples',sep=",")
sl = np.random.randint(len(glu_samples),size = 100)
gl = np.random.randint(len(gly_samples),size = 100)
params = ['k_4','KS_4','al_4','k_7','KS_7','al_7']
Z = np.swapaxes((sl,gl),0,1)
def sim_dif_params(Z):
    s,g = Z
    t_plot = np.linspace(df[paper][exp]['t_points'][0],df[paper][exp]['t_points'][-1],100)
    a = list(gly_samples.iloc[g]) + list(glu_samples.iloc[s])
    for p in params:
        pr[p] = a[params.index(p)]
    y,t_plot = nf.integrate_stiff(df[paper][exp]['y0'],t_plot,pr,I_f)
    return  y,t_plot


def plot_this():
    fig,ax = plt.subplots(1,1)
    ax.tick_params(labelsize=22)
    ax.set_ylim([0,1.5])
    ax.set_xlabel('Time ($d$)',fontsize=25)
    ax.set_ylabel('Concentration ($g.L^{-1}$)',fontsize=25)
    ax.tick_params(labelsize=25)

    plt.legend(fontsize=15)
    plt.show()

Ta_glu = 3.1 
Ta_gly = 1.4
Ta_gly_glu = 1.14
fr = 0.21
