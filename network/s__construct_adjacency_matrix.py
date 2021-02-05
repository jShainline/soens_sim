import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

from _functions_network import populate_hierarchy__power_law, populate_hierarchy__geometrical, generate_out_degree_distribution, generate_spatial_structure, determine_indices, neuron_level_rentian_scaling__with_spatial_dependence
from _plotting_network import plot_hierarchy__power_law, plot_hierarchy__geometrical, plot_out_degree_distribution, plot_node_degree_vs_space, plot_distance_matrix, plot_A

plt.close('all')

#%%

# this script produces the following data structures:
    # hierarchy
    # out_degree_distribution
    # spatial_information

#%% establish network hierarchy

print('constructing network hierarchy ... ')

# power-law hierarchy ( number of modules at level h of hierarchy: M_h = n_0 * h**(-gamma) )
num_levels_hier = 3
num_row_col_0 = 7
gamma = 2
# hierarchy = populate_hierarchy__power_law(num_row_col_0**2,num_levels_hier,float(gamma))
# plot_hierarchy__power_law(hierarchy)

# geometrical hierarchy
num_levels_hier = 4
num_row_col_0 = 7
delta_num_row_col = 2
hierarchy = populate_hierarchy__geometrical(num_row_col_0,delta_num_row_col,num_levels_hier)
plot_hierarchy__geometrical(hierarchy)

#%% generate out-degree distribution

print('generating out-degree distribution ... ')

degree_distribution__functional_form = 'gaussian' # 'gaussian' or 'power_law'
# note that we can let alpha = 0 in the power-law distribution for a random degree distribution or we can let st_dev = 0 in the gaussian distribution for a delta function (all nodes same degree)

if degree_distribution__functional_form == 'gaussian':
    
    # gaussian degree params
    center = 40 # mean of gaussian distribution
    st_dev = 5 # standard deviation of gaussian distribution    
    out_degree_distribution = generate_out_degree_distribution(degree_distribution__functional_form, center = center, st_dev = st_dev, num_nodes = hierarchy['total_num_nodes'])
    
elif degree_distribution__functional_form == 'power_law':
    
    # power-law degree params
    k_out_min = 2 # minimum degree allowed
    alpha = 2 # exponent, p(k) ~ k**(-alpha)
    out_degree_distribution = generate_out_degree_distribution(degree_distribution__functional_form, k_out_min = k_out_min, alpha = alpha, num_nodes = hierarchy['total_num_nodes'])
    
num_bins = 10    
plot_out_degree_distribution(out_degree_distribution,num_bins)

# avg_num_edges = 40

# total_num_nodes = hierarchy['total_num_nodes']
# tot_num_edges__target = avg_num_edges*total_nodes
# tot_num_edges__current = np.sum(out_degree_distribution['node_degrees'])

# diff = tot_num_edges__target-tot_num_edges__current
# print('diff = {}'.format(diff))
    

#%% place nodes in space

print('assigning spatial coordinates ... ')

spatial_information = generate_spatial_structure(hierarchy,out_degree_distribution)
    
plot_node_degree_vs_space(hierarchy,spatial_information)    
# plot_distance_matrix(spatial_information['distance_mat'],hierarchy['total_num_nodes'])
# plot_distance_matrix(spatial_information['distance_mat__corner'],hierarchy['total_num_nodes'])

#%% find indices of nodes within and external to each module at each level of hierarchy

# print('finding intra- and inter-modular indices ... ')

# indices_arrays, hierarchy = determine_indices(hierarchy,spatial_information)


#%% make connections

print('assigning edges based on neuron-level rent''s rule with spatial dependence ... ')

rentian_exponent = 0.8

spatial_information['spatial_dependence'] = 'exponential' # 'exponential' or 'power_law'
spatial_information['exponential_decay_factor'] = 3 # P_h(r_k = r) = A exp( -r/r_h ) where r_h = edf*nnrc_h where edf is the exponential_decay_factor set here and nnrc_h is hierarchy['num_nodes_row_col'][h], the length of a module in lattice constants at this level of hierarchy
# A, rent, out_degree_distribution = neuron_level_rentian_scaling__with_spatial_dependence(hierarchy,indices_arrays,spatial_information,out_degree_distribution,rentian_exponent)
A, rent, out_degree_distribution = neuron_level_rentian_scaling__with_spatial_dependence(hierarchy,spatial_information,out_degree_distribution,rentian_exponent)

plot_A(A)




#%% make connections

# exponential decay params
# decay_length = 5 # units of lattice constant

# num_bins = 20 # for plotting





#%% rescale out-degree distribution to achieve exact desired number of edges 

# avg_num_edges = 40

# total_nodes = hierarchy['total_nodes']
# tot_num_edges__target = avg_num_edges*total_nodes
# tot_num_edges__current = np.sum(out_degree_distribution['node_degrees'])

# diff = tot_num_edges__target-tot_num_edges__current
# print('diff = {}'.format(diff))

# scale_factor = tot_num_edges__current/tot_num_edges__target
# out_degree_distribution['node_degrees'] = np.round(out_degree_distribution['node_degrees']/scale_factor)
# tot_num_edges__current = np.sum(out_degree_distribution['node_degrees'])

# diff = tot_num_edges__target-tot_num_edges__current
# print('diff = {}'.format(diff))

# indices = np.random.randint(0,total_nodes-1,np.round(np.abs(diff)).astype(int))
# out_degree_distribution['node_degrees'][indices] += np.sign(diff) 
# tot_num_edges__current = np.sum(out_degree_distribution['node_degrees'])

# diff = tot_num_edges__target-tot_num_edges__current
# print('diff = {}'.format(diff))

# plot_out_degree_distribution(out_degree_distribution,num_bins)




