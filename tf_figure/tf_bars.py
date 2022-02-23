import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

tf_names = ['Batch','Batch GLU/GLY','Fedbatch GLU','Fedbatch GLY','Fedbatch GLU/GLY']
tf_list = [392/24,65/24,3.19,1.40,1.15]
tf_std = [2.2,0.39,0.35,0.21,0.19]
ind = np.arange(len(tf_names))
width = 0.35
# p1 = plt.bar(ind, tf_list, width, yerr=menStd)
fig, ax = plt.subplots()
p1 = ax.bar(ind, tf_list, width,yerr=tf_std)
ax.set_ylabel('$T_{f} \ (days)$',fontsize =30)
ax.set_xticks(ind,tf_names,fontsize = 30)
ax.tick_params(labelsize=25)
plt.show()
