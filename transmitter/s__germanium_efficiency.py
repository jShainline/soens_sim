#%%
import numpy as np
from matplotlib import pyplot as plt

from _functions import fermi_dirac
from _util import physical_constants
p = physical_constants()

plt.close('all')

#%% inputs

T = 4.2 # temperature in kelvin

E_f = 0 # fermi level
V = 0.01 # applied voltage in volts

E_L = E_f
E_gamma = E_f+0.01*p['eV']

fd_L = fermi_dirac(E_L,E_f,V,T)
fd_gamma = fermi_dirac(E_gamma,E_f,V,T)

efficiency = fd_gamma/(fd_L+fd_gamma)

print('Efficiency = {:6.4e}'.format(efficiency))

