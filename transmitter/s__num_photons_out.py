import numpy as np
from matplotlib import pyplot as plt
from _functions import Ljj

from _util import physical_constants, color_dictionary, set_plot_params
p = physical_constants()
colors = color_dictionary()

fig_size = set_plot_params('large') # 'publication' or 'large'

plt.close('all')

#%% inputs
eta = 1e-4 # LED total efficiency
tau_on = 50e-9 # time LED current is on (assuming square pulse)
I_LED = 10e-6 # LED bias current

print('With eta = {:6.4e}, tau_on = {:6.4f}ns, and I_LED = {:6.4f}uA, N_ph = {:6.4e}'.format(eta,tau_on*1e9,I_LED*1e6,eta*tau_on*I_LED/p['e']))