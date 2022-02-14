import numpy
import notebook_functions as nf

file_params = '../Parameters_inhibition.xlsx'
pr = nf.get_init_params(file_params)

v_M = {'BUT': 88,'ACE': 59 ,'X': 186, 'GAP':170, 'SUC': 118, 'GLU': 180, 'GLY': 92}
variables = ['BUT', 'ACE', 'SUC', 'X', 'GAP', 'GLU', 'GLY']