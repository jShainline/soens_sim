import numpy as np
from matplotlib import pyplot as plt

from _util import physical_constants, color_dictionary
p = physical_constants()
colors = color_dictionary()

plt.close('all')

#%% inputs

# assuming analog flux integration loops and a collection coil at every dendrite

N_vec = np.logspace(0,4,1000)

tn = 14 # width of plots in inches

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
# L_dc2 = 10e-12 # pH
L_dc3 = [100e-12,200e-12,400e-12] # pH

L_di2_vec = np.zeros([len(k1),len(L_dc3),len(N_vec)])

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
plt.suptitle('Max flux criterion\nIc = {:6.2f}uA, I_di_sat = {:6.2f}uA, L_dc1 = {:5.2f}pH, L_dc2 = L_dc3/10'.format(Ic*1e6,I_di_sat*1e6,L_dc1*1e12))
color_list = [['blue1','blue3','blue5'],['green1','green3','green5']]
linestyle_list = ['solid','dashdot']
for ii in range(len(k1)):
    for jj in range(len(L_dc3)):
        L_dc2 = L_dc3[jj]/10
        L_di2_vec[ii,jj,:] = (2*Ic/( k1[ii]**2 * k2[ii]**2 * L_dc1 * L_dc3[jj] * p['Phi0'])) * ( (p['Phi0']/(N_vec[:]*I_di_sat)) * (N_vec[:]*L_dc1 + L_dc2 + L_dc3[jj]) )**2
        ax.loglog(N_vec[:],L_di2_vec[ii,jj]*1e12, color = colors[color_list[ii][jj]], linestyle = linestyle_list[ii], label = 'k = {:4.2f}, L_dc3 = {:6.2f}pH'.format(k1[ii],L_dc3[jj]*1e12)) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])
# ax.loglog(N_vec,N_vec, linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  

ax.set_ylabel(r'$L_di2$ [pH]')
ax.set_xlabel(r'$N$')
ax.tick_params(axis = 'both')
ax.grid(which = 'both', axis = 'both')
ax.set_xlim([N_vec[0],N_vec[-1]])
ax.set_ylim([10,10000])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

#%% with inductances specified, calculate fraction of input activity required to drive dendrite to threshold

k1 = 0.5
k2 = 0.5

Idr_0_frac_vec = np.linspace(0,1,1000) # vector of fraction of JJ Ic at which squid is biased

P_over_N_vec = ((3*np.pi+2)/(2*np.pi)) *(1 - Idr_0_frac_vec)

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
plt.suptitle('Fraction of activity required for threshold\nIc = {:6.2f}uA, I_di_sat = {:6.2f}uA, k1 = {:4.2f}, k2 = {:4.2f}, L_dc1 = {:5.2f}pH, L_dc2 = L_dc3/10'.format(Ic*1e6,I_di_sat*1e6,k1,k2,L_dc1*1e12))
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

#%% dendritic tree - hierarchy analysis

N = 1000 # total number of synapses
n = 5 # fan-in factor (number of inputs to each dendrite)
h = np.log(N)/np.log(n)

h_vec = [3,4,5]
color_list = ['blue3','red3','green3']
fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
for ii in range(len(h_vec)):
    ax.loglog(N_vec,N_vec**(1/h_vec[ii]), color = colors[color_list[ii]], label = 'h = {:d}'.format(h_vec[ii]))
ax.set_ylabel(r'$n = N^{1/h}$')
ax.set_xlabel(r'$N$')
ax.tick_params(axis = 'both')
ax.grid(which = 'both', axis = 'both')
ax.set_xlim([N_vec[0],N_vec[-1]])
# ax.set_ylim([0.01,1])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

#%% dendritic tree

# make a function of h_vec

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
plt.suptitle('Num synapses with dendritic tree vs point neuron')
ax.plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**(h-1), color = colors['blue3']) # , label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj])

ax.set_ylabel(r'$P_dt/P_pn$')
ax.set_xlabel(r'$I_b/I_c$')
ax.tick_params(axis = 'both')
ax.grid(which = 'major', axis = 'both')
ax.set_xlim([0.4,1]) # Idr_0_frac_vec[0],Idr_0_frac_vec[-1]
ax.set_ylim([0,1])
# plt.subplots_adjust(wspace=0, hspace=0)

ax.legend()
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
plt.suptitle('Fraction of synapses required for threshold')
ax.plot(Idr_0_frac_vec[:],(P_over_N_vec[:])**h, color = colors['blue3'], label = 'DT')
ax.plot(Idr_0_frac_vec[:],(P_over_N_vec[:]), color = colors['red3'], label = 'PN')

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



#%% point-neuron fraction required to drive dendrite/neuron to threshold with collection coil


# P = np.zeros([len(Ic_vec),len(I0_th_frac_vec),len(Lsi2_vec),len(N_vec)])

# for ii in range(len(Ic_vec)):
#     for jj in range(len(I0_th_frac_vec)):
#         for kk in range(len(Lsi2_vec)):
#             I_si_sat = Ic_vec[ii]*I0_th_frac_vec[jj] # saturation current of SI loop
#             # I_si_sat = Ic_vec[ii] # saturation current of SI loop
#             P[ii,jj,kk,:] = ( (3*np.pi+2)/(np.sqrt(2)*np.pi) ) * ( np.sqrt(Ic_vec[ii]*p['Phi0'])/(k1*k2*I_si_sat) ) * (1-I0_th_frac_vec[jj]) * (N_vec*Lnc1 + Lnc2 + Lnc3_vec[kk]) / np.sqrt(Lnc1*Lnc3_vec[kk]*Lsi2_vec[kk])
    

# # color_list = ['blue3','red3','green3','yellow3']
# color_list = [['blue3','blue5'],['red3','red5'],['green3','green5']]
# linestyle_list = ['dashed','solid']
# for ii in range(len(Ic_vec)):
#     fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
#     plt.suptitle('Ic = {:5.2f}uA'.format(Ic_vec[ii]*1e6))
#     for jj in range(len(I0_th_frac_vec)):
#         for kk in range(len(Lsi2_vec)):
#             ax.loglog(N_vec,P[ii,jj,kk,:]/N_vec, linestyle = linestyle_list[kk], color = colors[color_list[jj][kk]], label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj]))
#     ax.loglog(N_vec,1/np.sqrt(N_vec), linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  
    
#     ax.set_ylabel(r'$P/N$')
#     ax.set_xlabel(r'$N$')
#     ax.tick_params(axis = 'both')
#     ax.grid(which = 'both', axis = 'both')
#     ax.set_xlim([N_vec[0],N_vec[-1]])
#     ax.set_ylim([0.01,1])
#     # plt.subplots_adjust(wspace=0, hspace=0)

#     ax.legend()
#     plt.tight_layout()
#     plt.show()

#%% cross talk

# active_fraction_vec = np.logspace([-3,0,1000])

# cross_talk_vec = np.zeros([len(L_b_vec),len(L_m_vec),len(N_vec)])

# for ii in range(len(L_b_vec)):
#     for jj in range(len(L_m_vec)):
#         epsilon = 2*L_m_vec[jj]/(N_vec[:]*L_s)
#         cross_talk_vec[ii,jj,:] = N_vec[:] * k**2 * L_m_vec[jj] / ( ( L_b_vec[ii]+L_m_vec[jj] ) * (1+epsilon[:]))

# fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))

# for ii in range(len(L_b_vec)):
#     for jj in range(len(L_m_vec)):
#         ax.loglog(N_vec,cross_talk_vec[ii,jj,:], linestyle = linestyle_list[jj], color = colors[color_list[ii][jj]], label = 'L_m = {:5.2f}pH, L_b = {:5.2f}nH'.format(L_m_vec[jj]*1e12,L_b_vec[ii]*1e9))
# # ax.loglog(N_vec,1/np.sqrt(N_vec), linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  

# ax.set_ylabel(r'$I^{si}_{ind}/I_0$')
# ax.set_xlabel(r'$N$')
# ax.tick_params(axis = 'both')
# ax.grid(which = 'both', axis = 'both')
# ax.set_xlim([N_vec[0],N_vec[-1]])
# # ax.set_ylim([0.01,1])
# # plt.subplots_adjust(wspace=0, hspace=0)

# ax.legend()
# plt.tight_layout()
# plt.show()

#%% dendritic tree



