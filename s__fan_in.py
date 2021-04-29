import numpy as np
from matplotlib import pyplot as plt
from _functions import Ljj

from _util import physical_constants, color_dictionary, set_plot_params
p = physical_constants()
colors = color_dictionary()

fig_size = set_plot_params('publication') # 'publication' or 'large'

plt.close('all')

#%% inputs

# assuming analog flux integration loops and a collection coil at every dendrite

N_vec = np.logspace(0,4,1000)

#%% choose Ic based on area-energy tradeoff 
Ic = 100e-6 # junction critical current

#%% design DR based on beta_L = 1
Ldr1_Ldr2 = p['Phi0']/(2*Ic)
Lj1 = p['Phi0']/(2*np.pi*Ic)
Lj2 = p['Phi0']/(4*Ic)

Ldr_tot = Ldr1_Ldr2+Lj1+Lj2
print('Ldr_tot = {:5.2f}pH = {:5.2f}pH'.format(Ldr_tot*1e12, 1e12*(p['Phi0']/Ic)*((3*np.pi+2)/(4*np.pi)) ))

#%% choose inductors based on max flux criterion
Phi_dr_max = p['Phi0']/2
I_di_sat = 0.4*Ic # Ic/2

k1 = [0.5,0.75]
k2 = [0.5,0.75]

L_dc1 = p['Phi0']/(2*Ic) # 10e-12 # pH
L_dc2_alpha = 0.1 # fraction of L_dc3 that is parasitic
L_dc3 = [100e-12,200e-12,400e-12] # pH

L_di2_vec = np.zeros([len(k1),len(L_dc3),len(N_vec)])

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (fig_size,0.75*fig_size))
plt.suptitle('Max flux criterion\nIc = {:6.2f}uA, I_di_sat = {:6.2f}uA, L_dc1 = {:5.2f}pH, L_dc2 = {:4.2f}L_dc3'.format(Ic*1e6,I_di_sat*1e6,L_dc1*1e12,L_dc2_alpha))
color_list = [['blue1','blue3','blue5'],['green1','green3','green5']]
linestyle_list = ['solid','dashdot']
for ii in range(len(k1)):
    for jj in range(len(L_dc3)):
        L_dc2 = L_dc2_alpha*L_dc3[jj]
        # L_di2_vec[ii,jj,:] = (2*Ic/( k1[ii]**2 * k2[ii]**2 * L_dc1 * L_dc3[jj] * p['Phi0'])) * ( (Phi_dr_max/(N_vec[:]*I_di_sat)) * (N_vec[:]*L_dc1 + L_dc2 + L_dc3[jj]) )**2
        L_di2_vec[ii,jj,:] = (1/L_dc3[jj]) * ( (Phi_dr_max/(k1[ii]*k2[ii]*N_vec[:]*I_di_sat)) * ( N_vec[:] + (2*Ic/p['Phi0'])*L_dc3[jj]*(1+L_dc2_alpha) ) )**2
        ax.loglog(N_vec[:],L_di2_vec[ii,jj,:]*1e12, color = colors[color_list[ii][jj]], linestyle = linestyle_list[ii], label = 'k = {:4.2f}, L_dc3 = {:6.2f}pH'.format(k1[ii],L_dc3[jj]*1e12)) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])
# ax.loglog(N_vec,N_vec, linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  

ax.set_ylabel(r'$L_di2$ [pH]')
ax.set_xlabel(r'$N$')
ax.tick_params(axis = 'both')
ax.grid(which = 'both', axis = 'both')
ax.set_xlim([N_vec[0],N_vec[-1]])
ax.set_ylim([1,10000])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

#%% with inductances specified, calculate fraction of input activity required to drive dendrite to threshold

k1 = 0.5
k2 = 0.5

Idr_0_frac_vec = np.linspace(0,1,1000) # vector of fraction of JJ Ic at which squid is biased

P_over_N_vec = ((3*np.pi+2)/(2*np.pi)) *(1 - Idr_0_frac_vec)

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
plt.suptitle('Fraction of activity required for threshold\nIc = {:6.2f}uA, I_di_sat = {:6.2f}uA, k1 = {:4.2f}, k2 = {:4.2f}, L_dc1 = {:5.2f}pH, L_dc2 = {:4.2f}L_dc3'.format(Ic*1e6,I_di_sat*1e6,k1,k2,L_dc1*1e12,L_dc2_alpha))
ax.plot(Idr_0_frac_vec[:],P_over_N_vec[:], color = colors['blue3']) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])

# ax.loglog(N_vec,np.sqrt(N_vec), linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  

ax.set_ylabel(r'$P/N$')
ax.set_xlabel(r'$I_b/I_c$')
ax.tick_params(axis = 'both')
ax.grid(which = 'major', axis = 'both')
ax.set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
ax.set_ylim([0,1])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()


#%% with inductances specified, calculate cross talk induced in one DI loop due to activity on others.

k1 = [0.5,0.75]
k2 = [0.5,0.75]

L_dc3 = 100e-12 # pH

L_di1 = [0,1e-9,10e-9]

I_ind_over_sat = np.zeros([len(k1),len(L_di1),len(N_vec)])

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
plt.suptitle('Cross talk when all loops are saturated\nIc = {:6.2f}uA, I_di_sat = {:6.2f}uA, L_dc1 = {:5.2f}pH, L_dc2 = {:4.2f}L_dc3, L_dc3 = {:6.2f}pH'.format(Ic*1e6,I_di_sat*1e6,L_dc1*1e12,L_dc2_alpha,L_dc3*1e12))
color_list = [['blue3','yellow3','green3'],['blue2','yellow2','green2']]
linestyle_list = ['solid','dotted']
for ii in range(len(k1)):
    for jj in range(len(L_di1)):
        L_di2 = (1/L_dc3) * ( (Phi_dr_max/(k1[ii]*k2[ii]*N_vec[:]*I_di_sat)) * ( N_vec[:] + (2*Ic/p['Phi0'])*L_dc3*(1+L_dc2_alpha) ) )**2
        L_di_tot = Ljj(Ic,0) + L_di1[jj] + L_di2
        L_dc_tot = N_vec[:]*L_dc1 + L_dc3*(1+L_dc2_alpha)
        # print('L_di2[0] = {:5.2f}pH, L_di2[-1] = {:5.2f}pH, L_di_tot[0] = {:5.2f}nH, L_di_tot[-1] = {:5.2f}nH'.format(L_di2[0]*1e12,L_di2[-1]*1e12,L_di_tot[0]*1e9,L_di_tot[-1]*1e9))
        I_ind_over_sat[ii,jj,:] = ( (k1[ii]**2 * L_di2 * L_dc1) / (L_di_tot*L_dc_tot) ) * N_vec[:]
        ax.loglog(N_vec[:],I_ind_over_sat[ii,jj,:], color = colors[color_list[ii][jj]], linestyle = linestyle_list[ii], label = 'k = {:4.2f}, L_di1 = {:5.2f}nH'.format(k1[ii],L_di1[jj]*1e9))

ax.set_ylabel(r'$I^{di}_{ind}/I^{di}_{sat}$')
ax.set_xlabel(r'$N$')
ax.tick_params(axis = 'both')
ax.grid(which = 'major', axis = 'both')
ax.set_xlim([N_vec[0],N_vec[-1]])
ax.set_ylim([1e-3,1])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

#%% dendritic tree - hierarchy analysis

N = 1000 # total number of synapses
n = 5 # fan-in factor (number of inputs to each dendrite)
h = np.log(N)/np.log(n)

h_vec = [2,3,4,5]
color_list = ['blue1','blue2','blue3','blue4','blue5']
fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
for ii in range(len(h_vec)):
    ax.loglog(N_vec,N_vec**(1/h_vec[ii]), color = colors[color_list[ii]], label = 'h = {:d}'.format(h_vec[ii]))
ax.set_ylabel(r'$n = N^{1/h}$')
ax.set_xlabel(r'$N$')
ax.tick_params(axis = 'both')
ax.grid(which = 'both', axis = 'both')
ax.set_xlim([N_vec[0],N_vec[-1]])
ax.set_ylim([1,1e2])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

#%% dendritic tree

# fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
# plt.suptitle('Num synapses with dendritic tree vs point neuron')

# color_list = ['blue3','red3','green3','yellow3']
# for ii in range(len(h_vec)):
#     ax.plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**(h_vec[ii]-1), color = colors[color_list[ii]], label = 'h = {:d}'.format(h_vec[ii])) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])

# ax.set_ylabel(r'$P_{dt}/P_{pn}$')
# ax.set_xlabel(r'$I_b/I_c$')
# ax.tick_params(axis = 'both')
# ax.grid(which = 'major', axis = 'both')
# ax.set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
# ax.set_ylim([0,1])
# # plt.subplots_adjust(wspace=0, hspace=0)

# ax.legend()
# plt.tight_layout()
# plt.show()

h_vec = [1,2,3,4,5]

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
plt.suptitle('Fraction of synapses required for threshold')
for ii in range(len(h_vec)):
    ax.plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**h_vec[ii], color = colors[color_list[ii]], label = 'DT, h = {:d}'.format(h_vec[ii]))

ax.set_ylabel(r'$P/N$')
ax.set_xlabel(r'$I_b/I_c$')
ax.tick_params(axis = 'both')
ax.grid(which = 'major', axis = 'both')
ax.set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
ax.set_ylim([0,1])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

# both plots

# fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = False, figsize = (fig_size,fig_size))
# plt.suptitle('Num synapses and fraction requried for threshold with dendritic tree vs point neuron')

# color_list = ['blue3','red3','green3','yellow3']
# for ii in range(len(h_vec)):
#     ax[0].plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**(h_vec[ii]-1), color = colors[color_list[ii]], label = 'h = {:d}'.format(h_vec[ii])) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])

# ax[0].set_ylabel(r'$P_{dt}/P_{pn}$')
# ax[0].tick_params(axis = 'both')
# ax[0].grid(which = 'major', axis = 'both')
# ax[0].set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
# ax[0].set_ylim([0,1])
# # plt.subplots_adjust(wspace=0, hspace=0)

# ax[0].legend()
# plt.tight_layout()
# plt.show()

# # fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False)
# # plt.suptitle('Fraction of synapses required for threshold')
# for ii in range(len(h_vec)):
#     ax[1].plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**h_vec[ii], color = colors[color_list[ii]], label = 'DT, h = {:d}'.format(h_vec[ii]))
# ax[1].plot(Idr_0_frac_vec[:],(P_over_N_vec[:]), color = colors['black'], label = 'PN')

# ax[1].set_ylabel(r'$P/N$')
# ax[1].set_xlabel(r'$I_b/I_c$')
# ax[1].tick_params(axis = 'both')
# ax[1].grid(which = 'major', axis = 'both')
# ax[1].set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
# ax[1].set_ylim([0,1])
# # plt.subplots_adjust(wspace=0, hspace=0)

# ax[1].legend()

# plt.subplots_adjust(wspace=0, hspace=0)
# # plt.tight_layout()
# plt.show()




