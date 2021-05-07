import numpy as np
from matplotlib import pyplot as plt
from _functions import Ljj

from _util import physical_constants, color_dictionary, set_plot_params
p = physical_constants()
colors = color_dictionary()

fig_size = set_plot_params('large') # 'publication' or 'large'

plt.close('all')

#%% inputs
Ic = 100e-6
L1 = 200e-12
L2 = 5.17e-12
k = 1
Ib = 140e-6
I_in = 10e-6

#%%
L_sq = 2*L2
L_sq_tot = 2*L2 + Ljj(Ic,86e-6) + Ljj(Ic,54e-6)
# L_sq_tot = 2*L2 # + Ljj(Ic,0) + Ljj(Ic,Ic)
I2__with_JJs = Ib/2 + (k*np.sqrt(L1*L2)/L_sq_tot)*I_in
I2__without_JJs = Ib/2 + (k*np.sqrt(L1*L2)/L_sq)*I_in

print('with JJs: I2 = {:7.4f}uA'.format(I2__with_JJs*1e6))
print('without JJs: I2 = {:7.4f}uA'.format(I2__without_JJs*1e6))