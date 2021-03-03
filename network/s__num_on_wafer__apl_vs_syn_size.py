import numpy as np
from matplotlib import pyplot as plt

from _util import physical_constants, color_dictionary
p = physical_constants()
colors = color_dictionary()

plt.close('all')

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14

plt.rcParams['figure.figsize'] = [14,11]
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 12

#%%
#--- inputs ---#
r_300 = 0.15 # radius of 300-mm wafer
w_wew = 1e-5 # pitch of wafer-edge waveguides
w_spd = 2.5e-5 # pitch of single-photon detectors in 1D
r_f125 = 125e-6/2 # radius of bare fiber
r_f250 = 250e-6/2 # radius of acrylic-clad fiber
gamma = 0.5772 # eulers constant

#%%
#--- calculated quantities ---#
h_8 = 2*r_300*np.sin(np.pi/8)
d_8 = 2*r_300*np.sqrt(1-np.sin(np.pi/8)**2)
A_8 = 2*np.sqrt(2)*r_300**2
A_4 = h_8**2
A_spd = w_spd**2
N_spd = A_8/A_spd
A_f125 = (3*np.sqrt(3)/2)*r_f125**2#for hexagonal packing
A_f250 = (3*np.sqrt(3)/2)*r_f250**2#for hexagonal packing
N_f125 = A_4/A_f125
N_f250 = A_4/A_f250

#%%
#--- output info ---#
print('N_wew = {:5.2e} wafer-edge waveguides on a 300-mm wafer cut into an octagon per cardinal direction'.format(h_8/w_wew))
print('N_wew = {:5.2e} total wafer-edge waveguides on a 300-mm wafer cut into an octagon with cardinal access'.format(4*h_8/w_wew))
print('N_spd = {:5.2e} total vertical photonic communication links on each face of a 300-mm wafer cut into an octagon'.format(N_spd))
print('N_f125 = {:5.2e} optical fibers if 125 um'.format(N_f125))
print('N_f250 = {:5.2e} optical fibers if 250 um'.format(N_f250))

#%% initial random graph as point of reference (Erdos-Renyi APL)

NVec = np.logspace(3,8,1000)
APL_vec = [2,3,4,5] # vector of average path lengths to consider
kExpectation_vec = [10,100,1000,10000]

APL_mat = np.zeros([len(NVec),len(kExpectation_vec)])
kExpectation_mat = np.zeros([len(NVec),len(APL_vec)])
    
for ii in range(len(kExpectation_vec)):
    APL_mat[:,ii] = ( np.log(NVec[:]) - gamma ) / np.log(kExpectation_vec[ii]) + 1/2
    
for ii in range(len(APL_vec)):
    kExpectation_mat[:,ii] = np.exp( ( np.log(NVec[:]) - gamma ) / ( APL_vec[ii] - 0.5 )  )


fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = False, figsize = (14,10))
fig.suptitle('Erdos-Renyi random graph')

color_list = ['blue3','red3','green3','yellow3']

for ii in range(len(kExpectation_vec)):
    ax[0].semilogx(NVec,APL_mat[:,ii], '-', color = colors[color_list[ii]], label = 'k = {:d}'.format(kExpectation_vec[ii]))
ax[0].set_ylabel(r'$\bar{L}$')
ax[0].legend()

for ii in range(len(APL_vec)):
    ax[1].loglog(NVec,kExpectation_mat[:,ii], '-', color = colors[color_list[ii]], label = 'APL = {:d}'.format(APL_vec[ii]))
ax[1].set_ylabel(r'$\bar{k}$')
ax[1].legend()

ax[1].set_xlabel(r'Total Nodes, $N_{tot}$')

ax[0].grid(which = 'both', axis = 'both')
ax[1].grid(which = 'both', axis = 'both')
ax[1].set_xlim([NVec[0],NVec[-1]])
ax[1].set_ylim([1,2e5])
plt.show()
plt.subplots_adjust(wspace=0, hspace=0)

#%% average path length versus number of neurons on wafer
n_planes__photonic = 20
n_planes__electronic = 10

w_wg__list = [1.5e-6,3e-6,6e-6]
w_sy__list = [10e-6,30e-6,90e-6]

N_300__vec = np.logspace(5,7,100)

L_apl__mat__photonic = np.zeros([len(N_300__vec),len(w_wg__list)])
for ii in range(len(w_wg__list)):
    
    A_n = A_8/N_300__vec
    
    k_bar__photonic = n_planes__photonic*np.sqrt(A_n)/w_wg__list[ii]
    L_apl__mat__photonic[:,ii] = ( np.log(N_300__vec) - gamma ) / ( np.log(k_bar__photonic) ) + 0.5
 
L_apl__mat__electronic = np.zeros([len(N_300__vec),len(w_sy__list)])
for ii in range(len(w_sy__list)):
    
    A_n = A_8/N_300__vec
    
    k_bar__electronic = n_planes__electronic*A_n/(w_sy__list[ii]**2)
    L_apl__mat__electronic[:,ii] = ( np.log(N_300__vec) - gamma ) / ( np.log(k_bar__electronic) ) + 0.5

color_list = ['blue3','green3','yellow3']
fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = False, figsize = (14,10))
plt.suptitle('Average path length versus number of neurons on 300-mm wafer\nnum_planes_photonic = {:d}; num_planes_electronic = {:d}'.format(n_planes__photonic,n_planes__electronic))

for ii in range(len(w_wg__list)):
    ax[0].semilogx(N_300__vec,L_apl__mat__photonic[:,ii], '-', color = colors[color_list[ii]], label = 'w_wg = {:3.1f}um'.format(w_wg__list[ii]*1e6))     

for ii in range(len(w_wg__list)):
    ax[1].semilogx(N_300__vec,L_apl__mat__electronic[:,ii], '-', color = colors[color_list[ii]], label = 'w_sy = {:3.1f}um'.format(w_sy__list[ii]*1e6))

ax[0].set_ylabel(r'Average path length, $\bar{L}$ (photonics-limited)') 
ax[0].legend() 

ax[1].set_xlabel(r'Num neurons on 300-mm wafer, $N_{300}$')
ax[1].set_ylabel(r'Average path length, $\bar{L}$ (electronics-limited)') 
ax[1].set_xlim([N_300__vec[0],N_300__vec[-1]])
ax[1].legend()
    
ax[0].grid(which = 'both', axis = 'both')
ax[1].grid(which = 'both', axis = 'both')
 

plt.subplots_adjust(wspace=0, hspace=0)
plt.show() 


#%% out degree versus number of neurons on wafer
n_planes__photonic = 20
n_planes__electronic = 10

w_wg__list = [1.5e-6,3e-6,6e-6]
w_sy__list = [10e-6,30e-6,90e-6]

N_300__vec = np.logspace(5,7,100)

k_bar__mat__photonic = np.zeros([len(N_300__vec),len(w_wg__list)])
for ii in range(len(w_wg__list)):
    
    A_n = A_8/N_300__vec
    
    k_bar__mat__photonic[:,ii] = n_planes__photonic*np.sqrt(A_n)/w_wg__list[ii]
 
k_bar__mat__electronic = np.zeros([len(N_300__vec),len(w_sy__list)])
for ii in range(len(w_sy__list)):
    
    A_n = A_8/N_300__vec
    
    k_bar__mat__electronic[:,ii] = n_planes__electronic*A_n/(w_sy__list[ii]**2)

color_list = ['blue3','green3','yellow3']
fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = False, figsize = (14,10))
plt.suptitle('Out degree versus number of neurons on 300-mm wafer\nnum_planes_photonic = {:d}; num_planes_electronic = {:d}'.format(n_planes__photonic,n_planes__electronic))

for ii in range(len(w_wg__list)):
    ax[0].loglog(N_300__vec,k_bar__mat__photonic[:,ii], '-', color = colors[color_list[ii]], label = 'w_wg = {:3.1f}um'.format(w_wg__list[ii]*1e6))     

for ii in range(len(w_wg__list)):
    ax[1].loglog(N_300__vec,k_bar__mat__electronic[:,ii], '-', color = colors[color_list[ii]], label = 'w_sy = {:3.1f}um'.format(w_sy__list[ii]*1e6))

ax[0].set_ylabel(r'Out degree, $\bar{k}$ (photonics-limited)') 
ax[0].legend() 

ax[1].set_xlabel(r'Num neurons on 300-mm wafer, $N_{300}$')
ax[1].set_ylabel(r'Out degree, $\bar{k}$ (electronics-limited)') 
ax[1].set_xlim([N_300__vec[0],N_300__vec[-1]])
ax[1].legend()
    
ax[0].grid(which = 'both', axis = 'both')
ax[1].grid(which = 'both', axis = 'both')
 

plt.subplots_adjust(wspace=0, hspace=0)
plt.show() 
 

#%% average path length vs synapse width

N_n__list = [1e5,1e6,1e7]
w_syn__vec = np.logspace(-6,-3,100)
A_syn = w_syn__vec**2
N_syn = A_8/A_syn

L_bar = np.zeros([len(w_syn__vec),len(N_n__list)])
indices_list = []
L_bar__10 = np.zeros([len(w_syn__vec),len(N_n__list)])
indices_list__10 = []

for ii in range(len(N_n__list)):
    # L_bar[:,ii] = ( np.log(N_n__list[ii]) - gamma ) / ( np.log( A_8/(w_syn__vec**2*N_n__list[ii]) ) ) + 1/2
    L_bar[:,ii] = ( np.log(N_n__list[ii]) - gamma ) / ( np.log(A_8) - np.log(N_n__list[ii]) - 2*np.log(w_syn__vec) ) + 1/2
    indices_list.append( np.where( L_bar[:,ii] > 0 ) )
    L_bar__10[:,ii] = ( np.log(N_n__list[ii]) - gamma ) / ( np.log(A_8) - np.log(N_n__list[ii]) - 2*np.log(w_syn__vec/np.sqrt(10)) ) + 1/2
    indices_list__10.append( np.where( L_bar__10[:,ii] > 0 ) )


color_list = ['blue1','green1','yellow1']
color_list__10 = ['blue3','green3','yellow3']
fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True, sharey = False, figsize = (14,10))
# plt.suptitle('Network average path length versus width of synapse')

for ii in range(len(N_n__list)):
    ax[0].semilogx(w_syn__vec[indices_list[ii]]*1e6,L_bar[indices_list[ii][0],ii], '-', color = colors[color_list[ii]], label = 'N_n = {:4.1e}, 1 plane'.format(N_n__list[ii]))    
    ax[0].semilogx(w_syn__vec[indices_list__10[ii]]*1e6,L_bar__10[indices_list__10[ii][0],ii], '-', color = colors[color_list__10[ii]], label = 'N_n = {:4.1e}, 10 planes'.format(N_n__list[ii]))     

ax[0].set_xlim([w_syn__vec[0]*1e6,w_syn__vec[-1]*1e6])
ax[0].set_ylim([1,3])
ax[0].set_ylabel(r'Network average path length, $\bar{L}$') 
ax[0].grid(which = 'both', axis = 'both')
ax[0].legend() 

#%% ratio of waveguide area to electronic synapse circuit area

p_e = 10 # number of planes of active electronic devices
p_p = 20 # number of planes of photonic routing waveguides
L_apl = 2 # average network path length
w_wg = 1.5e-6
w_sy__vec = np.logspace(-6,-3,1000)

kappa__mat = np.zeros([len(w_sy__vec),len(N_n__list)])
for ii in range(len(N_n__list)):
    kappa__mat[:,ii] = (p_e/p_p**2) * (w_wg/w_sy__vec[:])**2 * np.exp( ( np.log(N_n__list[ii]) - gamma ) / (L_apl-1/2) )
 
# fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (14,10))
# plt.suptitle('Ratio of waveguide area to electronic synapse circuit area\nw_wg = {:3.1f}um; L_apl = {:d}'.format(w_wg*1e6,L_apl))
color_list = ['blue3','green3','yellow3']
for ii in range(len(N_n__list)):
    ax[1].loglog(w_sy__vec*1e6,kappa__mat[:,ii], '-', color = colors[color_list[ii]], label = 'N_n = {:4.1e}, p_e = {:d}; p_p = {:d}; L_apl = {:d}'.format(N_n__list[ii],p_e,p_p,L_apl))    

ax[1].grid('on', which = 'both', axis = 'both')
ax[1].set_xlim([1,100])
ax[1].set_ylim([1e-1,1e2])
ax[1].set_xlabel(r'Width of electronic synapse [$\mu$m]')
ax[1].set_ylabel(r'Ratio of waveguide area to electronic synapse circuit area, $\kappa$') 
ax[1].legend()  

plt.subplots_adjust(wspace=0, hspace=0)
 
plt.show()   


#%% scraps

# connection fraction
# N_n__list = [1e5,1e6,1e7]
# w_syn__vec = np.logspace(-6,-3,100)
# A_syn = w_syn__vec**2
# N_syn = A_8/A_syn
# fig, ax = plt.subplots(1,1, figsize = (14,10))
# plt.suptitle('fraction of possible synaptic connections versus width of synapse')
# color_list = ['blue1','red1','green1','yellow1']
# for ii in range(len(N_n__list)):
#     k_out = N_syn/N_n__list[ii]
#     k_frac = N_n__list[ii]*k_out/(N_n__list[ii]**2) # fraction of possible synaptic connections
#     ax.loglog(w_syn__vec*1e6,k_frac, '-', color = colors[color_list[ii]], label = 'N_n = {:4.2e}, 1 plane'.format(N_n__list[ii]))  
# color_list = ['blue3','green3','yellow3']
# for ii in range(len(N_n__list)):
#     k_out = 10*N_syn/N_n__list[ii]
#     k_frac = N_n__list[ii]*k_out/(N_n__list[ii]**2) # fraction of possible synaptic connections
#     ax.loglog(w_syn__vec*1e6,k_frac, '-', color = colors[color_list[ii]], label = 'N_n = {:4.2e}, 10 planes'.format(N_n__list[ii]))     
# ax.grid('on')
# ax.set_xlim([w_syn__vec[0]*1e6,w_syn__vec[-1]*1e6])
# ax.set_ylim([1/N_n__list[-1],0.1])
# ax.set_xlabel(r'Width of synapse [$\mu$m]')
# ax.set_ylabel(r'Connection fraction') 
# ax.legend()  
# plt.show()  

