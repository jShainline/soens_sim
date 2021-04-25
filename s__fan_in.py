import numpy as np
from matplotlib import pyplot as plt

from _util import physical_constants, color_dictionary
p = physical_constants()
colors = color_dictionary()

plt.close('all')

#%% inputs

Ic_vec = np.asarray([40e-6,70e-6,100e-6]) # vector of JJ Ics
I0_th_frac_vec = np.asarray([0.7,0.9]) # vector of fraction of JJ Ic at which squid is biased

Lsi1 = 100e-9 # self-inductance of SI loop
Lsi2_vec = np.asarray([200e-12,400e-12]) # output inductance from synapse/dendrite to collection coil

Lnc1 = 10e-12 # input inductance to collection coil
Lnc2 = 100e-12 # parasitic inductance of collection coil
Lnc3_vec = Lsi2_vec # output inductance from collection to squid of next stage

k1 = 0.75 # MI coupling factor between synapses/dendrites and collection coil
k2 = 0.75 # MI coupling factor between collection coil and squid

N_vec = np.logspace(0,4,1000)

tn = 14 # width of plots in inches


#%% point-neuron fraction required to drive dendrite/neuron to threshold with collection coil


P = np.zeros([len(Ic_vec),len(I0_th_frac_vec),len(Lsi2_vec),len(N_vec)])

for ii in range(len(Ic_vec)):
    for jj in range(len(I0_th_frac_vec)):
        for kk in range(len(Lsi2_vec)):
            I_si_sat = Ic_vec[ii]*I0_th_frac_vec[jj] # saturation current of SI loop
            # I_si_sat = Ic_vec[ii] # saturation current of SI loop
            P[ii,jj,kk,:] = ( (3*np.pi+2)/(np.sqrt(2)*np.pi) ) * ( np.sqrt(Ic_vec[ii]*p['Phi0'])/(k1*k2*I_si_sat) ) * (1-I0_th_frac_vec[jj]) * (N_vec*Lnc1 + Lnc2 + Lnc3_vec[kk]) / np.sqrt(Lnc1*Lnc3_vec[kk]*Lsi2_vec[kk])
    

# color_list = ['blue3','red3','green3','yellow3']
color_list = [['blue3','blue5'],['red3','red5'],['green3','green5']]
linestyle_list = ['dashed','solid']
for ii in range(len(Ic_vec)):
    fig, ax = plt.subplots(nrows = 1, ncols = 1, sharex = True, sharey = False, figsize = (tn,tn/1.618))
    plt.suptitle('Ic = {:5.2f}uA'.format(Ic_vec[ii]*1e6))
    for jj in range(len(I0_th_frac_vec)):
        for kk in range(len(Lsi2_vec)):
            ax.loglog(N_vec,P[ii,jj,kk,:]/N_vec, linestyle = linestyle_list[kk], color = colors[color_list[jj][kk]], label = 'Lsi2 = Lnc3 = {:5.2f}pH, I0_th_frac_vec = {:5.2f}'.format(Lsi2_vec[kk]*1e12,I0_th_frac_vec[jj]))
    ax.loglog(N_vec,1/np.sqrt(N_vec), linestyle = 'dotted', color = colors['grey9'], label = '1/sqrt(N)')  
    
    ax.set_ylabel(r'$P/N$')
    ax.set_xlabel(r'$N$')
    ax.tick_params(axis = 'both')
    ax.grid(which = 'both', axis = 'both')
    ax.set_xlim([N_vec[0],N_vec[-1]])
    ax.set_ylim([0.01,1])
    # plt.subplots_adjust(wspace=0, hspace=0)

    ax.legend()
    plt.tight_layout()
    plt.show()

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

N = 1000 # total number of synapses
n = 10 # fan-in factor (number of inputs to each dendrite)
h = np.log(N)/np.log(n)