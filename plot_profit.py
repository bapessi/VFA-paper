import matplotlib.pyplot as plt
import max_profit_onlygly as mp
import max_profit_onlyglu as mglu
import numpy as np
ta_list = np.arange(0,24*10,0.1)

profit_list = []
for ta in ta_list:
    # profit = mp.get_profit_fedbatch(ta)
    profit = mglu.get_profit_fedbatch(ta)
    profit_list.append(profit)

plt.plot(ta_list,profit_list)
plt.show()
