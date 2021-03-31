import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import copy
import time
# from matplotlib import pyplot as plt
import powerlaw
import networkx as nx
import igraph as ig
# from networkx.algorithms import community
import community as community_louvain

from _util import color_dictionary # physical_constants
colors = color_dictionary()

#%%

def A_random(num_nodes,num_edges):
    
    A = np.zeros([num_nodes,num_nodes])
    connections = np.random.randint(0,num_nodes,size=[2,num_edges])
    print('shape(connections) = {}'.format(np.shape(connections)))
    print(connections)
    for ii in range(num_edges):
        A[connections[0,ii],connections[1,ii]] = 1
    # print('shape(A) = {}'.format(np.shape(A)))    
    
    return A

def populate_hierarchy__power_law(num_nodes_0,num_levels_hier,gamma, plot = True): # power-law hierarchy ( number of modules at level h of hierarchy: M_h = n_0 * h**(-gamma) )
    
    num_modules_list = np.zeros([num_levels_hier])
    n_h__num_nodes_vs_hierarchy = np.zeros([num_levels_hier])
    n_h__num_nodes_vs_hierarchy[0] = 1
    
    h_vec = np.arange(1,num_levels_hier+1,1)
    for h in h_vec:
        num_modules_list[h-1] = np.round( num_nodes_0 * h**(-gamma) ) # num_nodes_0 * h**(-gamma) # 
        n_h__num_nodes_vs_hierarchy[h-1] = ( num_nodes_0**h )*( np.prod(np.arange(1,h+1,1)**(-gamma)) ) # num_nodes_0*( np.prod(num_nodes_0*np.arange(2,h+1,1)**(-gamma)) ) # 
    
    num_nodes_per_module = n_h__num_nodes_vs_hierarchy/num_modules_list    
    total_nodes = n_h__num_nodes_vs_hierarchy[-1]
    inter_modular_nodes = n_h__num_nodes_vs_hierarchy*(1-1/num_modules_list)

    nnrc = []
    nmrc = []
    for ii in range(h['num_levels_hier']):
        nnrc.append( np.sqrt(num_nodes_per_module[ii]) ) # num nodes in each row/col within a module at the lower level of hierarchy
        nmrc.append( np.sqrt(num_modules_list[ii]) ) # num modules in each row/col at this level of hierarchy
        
    hierarchy = dict()
        
    hierarchy['num_nodes_0'] = num_nodes_0
    hierarchy['num_levels_hier'] = num_levels_hier
    hierarchy['gamma'] = gamma
    hierarchy['h_vec'] = h_vec
    
    hierarchy['num_modules_list'] = num_modules_list # number of modules at each level of hierarchy
    hierarchy['n_h__num_nodes_vs_hierarchy'] = n_h__num_nodes_vs_hierarchy # number of nodes at each level of hierarchy
    hierarchy['num_nodes_per_module'] = num_nodes_per_module # number of nodes per module at each level of hierarchy
    hierarchy['inter_modular_nodes'] = inter_modular_nodes # number of inter-modular nodes at each level of hierarchy
    hierarchy['total_num_nodes'] = total_nodes # total number of nodes in the network 
    hierarchy['num_nodes_row_col'] = nnrc
    hierarchy['num_modules_row_col'] = nmrc        

    return hierarchy

def populate_hierarchy__geometrical(num_row_col_0 = 9, delta_num_row_col = 2, num_levels_hier = 4, plot = True): # power-law hierarchy ( number of modules at level h of hierarchy: M_h = n_0 * h**(-gamma) )
    
    H = num_levels_hier
    M_h__num_submodules_vs_hierarchy = np.zeros([H])
    n_h__num_nodes_vs_hierarchy = np.zeros([H])
    # n_h__num_nodes_vs_hierarchy[0] = 1
    
    num_nodes_0 = num_row_col_0**2
    a = delta_num_row_col
    b = num_row_col_0+delta_num_row_col
    
    h_vec = np.arange(1,H+1,1)
    # print('h_vec = {}'.format(h_vec))
    for h in h_vec:
        # M_h__num_modules_vs_hierarchy[h-1] = np.round( (11 - 2*h)**2 )
        M_h__num_submodules_vs_hierarchy[h-1] = np.round( ( b-a*h )**2 )
        # print('np.prod( M_h__num_modules_vs_hierarchy[1:h] ) = {}'.format(np.prod( M_h__num_modules_vs_hierarchy[1:h] )))
        n_h__num_nodes_vs_hierarchy[h-1] = num_nodes_0*np.prod( M_h__num_submodules_vs_hierarchy[1:h] )
    
    num_nodes_per_module = n_h__num_nodes_vs_hierarchy/M_h__num_submodules_vs_hierarchy    
    total_nodes = n_h__num_nodes_vs_hierarchy[-1]
    inter_modular_nodes = n_h__num_nodes_vs_hierarchy*(1-1/M_h__num_submodules_vs_hierarchy)
        
    nnrc = np.zeros(H,int)
    nmrc = np.zeros(H,int)
    mh = np.zeros(H,int)
    for h in range(H):
        nnrc[h] = np.sqrt(num_nodes_per_module[h]) # num nodes in each row/col within a module at the lower level of hierarchy
        nmrc[h] = np.sqrt(M_h__num_submodules_vs_hierarchy[h]) # num modules in each row/col at this level of hierarchy
        mh[h] = np.prod(M_h__num_submodules_vs_hierarchy[h:H-1])
        
    hierarchy = dict()
    
    hierarchy['num_nodes_0'] = num_nodes_0
    hierarchy['H__num_levels_hier'] = num_levels_hier
    hierarchy['h_vec'] = h_vec
    
    hierarchy['M_h__num_submodules_vs_hierarchy'] = M_h__num_submodules_vs_hierarchy.astype(int) # number of sub-modules within each module at each level of hierarchy
    hierarchy['m_h__num_modules_vs_hierarchy'] = mh # number of modules at each level of hierarchy    
    # hierarchy['num_sub_modules_list'] = num_sub_modules_list.astype(int) # number of sub-modules contained in a module at each level of hierarchy above the lowest
    # hierarchy['map_to_upper'] = map_to_upper.astype(int) # index of containing module next level up in hierachy for each level of hierarchy except the highest
    hierarchy['n_h__num_nodes_vs_hierarchy'] = n_h__num_nodes_vs_hierarchy.astype(int) # number of nodes at each level of hierarchy
    hierarchy['num_nodes_per_module'] = num_nodes_per_module.astype(int) # number of nodes per module at each level of hierarchy
    hierarchy['inter_modular_nodes'] = inter_modular_nodes.astype(int) # number of inter-modular nodes at each level of hierarchy
    hierarchy['N_h__total_num_nodes'] = total_nodes.astype(int) # total number of nodes in the network         
    # hierarchy['num_row_col'] = np.sqrt(total_nodes).astype(int) # number of rows and columns in the square network    
    hierarchy['num_nodes_row_col'] = nnrc
    hierarchy['num_modules_row_col'] = nmrc  

    return hierarchy

def populate_hierarchy__arbitrary_square(num_sub_modules__row_col__list = [4,4,4]): # power-law hierarchy ( number of modules at level h of hierarchy: M_h = n_0 * h**(-gamma) )
    
    hierarchy = dict()
    
    H = len(num_sub_modules__row_col__list)
    M_h__num_submodules_vs_hierarchy = (np.asarray(num_sub_modules__row_col__list)**2).astype(int)
    n_h__num_nodes_per_module_vs_hierarchy = np.zeros(H)
    for h in range(H):
        n_h__num_nodes_per_module_vs_hierarchy[h] = np.prod(M_h__num_submodules_vs_hierarchy[0:h])
        
    num_nodes_per_module = n_h__num_nodes_per_module_vs_hierarchy
    total_nodes = n_h__num_nodes_per_module_vs_hierarchy[-1]
        
    nnrc = np.zeros(H,int)
    nmrc = np.zeros(H,int)
    mh = np.zeros(H,int)
    for h in range(H):
        nnrc[h] = np.sqrt(num_nodes_per_module[h]) # num nodes in each row/col within a module at the lower level of hierarchy
        nmrc[h] = np.sqrt(M_h__num_submodules_vs_hierarchy[h]) # num modules in each row/col at this level of hierarchy
        mh[h] = np.prod(M_h__num_submodules_vs_hierarchy[h:H-1])
        
    hierarchy['H__num_levels_hier'] = H
    hierarchy['h_vec'] = np.arange(1,H+1,1)
    hierarchy['M_h__num_submodules_vs_hierarchy'] = M_h__num_submodules_vs_hierarchy # number of sub-modules within each module at each level of hierarchy
    hierarchy['m_h__num_modules_vs_hierarchy'] = mh # number of modules at each level of hierarchy    
    hierarchy['n_h__num_nodes_per_module_vs_hierarchy'] = n_h__num_nodes_per_module_vs_hierarchy.astype(int) # number of nodes at each level of hierarchy
    hierarchy['N_h__total_num_nodes'] = total_nodes.astype(int) # total number of nodes in the network            
    hierarchy['num_nodes_row_col'] = nnrc
    hierarchy['num_modules_row_col'] = nmrc  

    return hierarchy

def generate_out_degree_distribution(out_degree_functional_form = 'power-law', **kwargs):
    
    if 'num_nodes' in kwargs:
        num_nodes = kwargs['num_nodes']
    else:
        raise ValueError('[_functions_network/generate_degree_distribution] You must specify the number of nodes in the network (num_nodes)') 
    
    out_degree_distribution = dict()
    if out_degree_functional_form == 'gaussian':

        if 'center' in kwargs:
            center = kwargs['center']
        else: 
            raise ValueError('[_functions_network/generate_degree_distribution] For a gaussian out-degree distribution, you must specify the mean of the gaussian distribution (center)')
        if 'st_dev' in kwargs:
            st_dev = kwargs['st_dev']
        else: 
            raise ValueError('[_functions_network/generate_degree_distribution] For a gaussian out-degree distribution, you must specify the standard deviation (st_dev)')        
        
        out_degree_distribution['functional_form'] = 'gaussian'
        out_degree_distribution['center'] = center
        out_degree_distribution['st_dev'] = st_dev
        out_degree_distribution['num_nodes'] = int(num_nodes)
        out_degree_distribution['node_degrees'] = np.flipud(np.round(np.sort(np.random.normal(center,st_dev,int(num_nodes))))) # gaussian degree distribution
        
    if out_degree_functional_form == 'power_law':
        
        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
        else:
            raise ValueError('[_functions_network/generate_degree_distribution] For a power-law out-degree distribution, you must specify the exponent (alpha)')
        if 'k_out_min' in kwargs:
            k_out_min = kwargs['k_out_min']
        else:
            raise ValueError('[_functions_network/generate_degree_distribution] For a power-law out-degree distribution, you must specify the minimum out degree (k_out_min)')           
        
        out_degree_distribution['functional_form'] = 'power_law'
        out_degree_distribution['k_out_min'] = k_out_min
        out_degree_distribution['alpha'] = alpha
        out_degree_distribution['num_nodes'] = int(num_nodes)
        out_degree_distribution['node_degrees'] = np.flipud(np.round(np.sort(powerlaw.Power_Law(xmin = k_out_min, parameters = [alpha]).generate_random(int(num_nodes))))) # power-law degree distribution
        
    return out_degree_distribution

def generate_spatial_structure(hier,o_d_d):

    spatial_information = dict()
        
    H = hier['H__num_levels_hier']
    Nh = hier['N_h__total_num_nodes']
    nh = hier['n_h__num_nodes_per_module_vs_hierarchy']
    Mh = hier['M_h__num_submodules_vs_hierarchy']
    mh = hier['m_h__num_modules_vs_hierarchy']
    nmrc = hier['num_modules_row_col']
    nnrc = hier['num_nodes_row_col']
    nnpm = nh # hier['num_nodes_per_module']

    # first establish modules and available sites for nodes within modules, all corner referenced to bottom left   
    mod_center_coords = []
    mod_corner_coords = []
    mod_coords_x = []
    mod_coords_y = []
    h_increment = 0.1
    for h in range(H):
        mod_center_coords.append([])
        mod_corner_coords.append([])
        mod_coords_x.append([])
        mod_coords_y.append([])
        for m in range(mh[h]):
            x_mod = np.floor(m/np.prod(nmrc[h:H-1]))
            y_mod = m-np.prod(nmrc[h:H-1])*x_mod
            delta = np.round( (nnrc[h]-0.1)/2 )
            mod_center_coords[h].append([x_mod*nnrc[h]+delta,y_mod*nnrc[h]+delta])
            mod_corner_coords[h].append([x_mod*nnrc[h],y_mod*nnrc[h]])
            mod_coords_x[h].append([x_mod*nnrc[h]-h*h_increment,x_mod*nnrc[h]+nnrc[h]-1+h*h_increment,x_mod*nnrc[h]+nnrc[h]-1+h*h_increment,x_mod*nnrc[h]-h*h_increment,x_mod*nnrc[h]-h*h_increment])
            mod_coords_y[h].append([y_mod*nnrc[h]-h*h_increment,y_mod*nnrc[h]-h*h_increment,y_mod*nnrc[h]+nnrc[h]-1+h*h_increment,y_mod*nnrc[h]+nnrc[h]-1+h*h_increment,y_mod*nnrc[h]-h*h_increment])
    spatial_information['mod_center_coords'] = mod_center_coords 
    spatial_information['mod_corner_coords'] = mod_corner_coords 
    spatial_information['mod_coords_x'] = mod_coords_x 
    spatial_information['mod_coords_y'] = mod_coords_y 
    
    # find order of modules descending from closest to center at each level of hierarchy
    module_order = []
    center = np.sqrt(Nh)/2    
    for h in range(H):
        module_order.append([])
        distance_list = np.zeros(mh[h])
        for m in range(mh[h]):
            distance_list[m] = ( (mod_center_coords[h][m][0] - center)**2 + (mod_center_coords[h][m][1]-center)**2 )**(1/2) # euclidean distance
        module_order[h] = np.argsort(distance_list)
    spatial_information['module_order'] = module_order
    
    # determine intra-module indices and distance from center
    intra_modular_indices = []   
    distance_from_center = []
    for h in range(H):
        intra_modular_indices.append([])
        distance_from_center.append([])
        for m in range(mh[h]):
            intra_modular_indices[h].append([])
            distance_from_center[h].append([])
            temp_vec1 = np.zeros(nnpm[h])
            temp_vec2 = np.zeros(nnpm[h])
            for nn in range(nnpm[h]):
                dx = np.floor(nn/nnrc[h])
                dy = nn - nnrc[h]*dx
                x = mod_corner_coords[h][m][0] + dx
                y = mod_corner_coords[h][m][1] + dy
                temp_vec1[nn] = nnrc[-1]*x+y
                temp_vec2[nn] = ( ( mod_center_coords[h][m][0] - x )**2 + ( mod_center_coords[h][m][1] - y )**2 )**(1/2)
            intra_modular_indices[h][m] = temp_vec1
            distance_from_center[h][m] = temp_vec2
    spatial_information['intra_modular_indices'] = intra_modular_indices                
    spatial_information['distance_from_center'] = distance_from_center                
                
    # step through modules, adding one node to unoccupied site closest to center
    available_sites__distances = copy.deepcopy(distance_from_center)
    available_sites__indices = copy.deepcopy(intra_modular_indices)
    h_vec = np.flipud(np.arange(1,H,1))
    node_coords = []
    node_indices = []
    num_nodes_placed = 0
    while num_nodes_placed <= Nh:
        for h in h_vec:
            for m in module_order[h]:
                if num_nodes_placed < Nh:
                    mindex = np.asarray(available_sites__distances[h][m]).argmin()
                    ii = available_sites__indices[h][m][mindex]
                    if ii not in node_indices:
                        node_indices.append(ii)
                        x = np.floor(ii/nnrc[-1])
                        node_coords.append(np.array([x,ii-nnrc[-1]*x]).astype(int))
                        num_nodes_placed += 1
                    available_sites__distances[h][m] = np.delete(available_sites__distances[h][m],mindex)
                    available_sites__indices[h][m] = np.delete(available_sites__indices[h][m],mindex)
                else:
                    num_nodes_placed += 1
    spatial_information['node_coords'] = node_coords                
    spatial_information['node_indices'] = node_indices                
    
    # generate out_degree vs spatial coords
    degree_xy = np.zeros([nnrc[-1],nnrc[-1]])
    for ii in range(Nh):
        degree_xy[node_coords[ii][0],node_coords[ii][1]] = o_d_d['node_degrees'][ii] # o_d_d['node_degrees'][node_indices[ii].astype(int)]                
    spatial_information['degree_xy'] = degree_xy   
           
    # generate map from each module to the containing module in the next level of hierarcy
    map_to_upper = []
    for ii_sub in range(Nh):
        map_to_upper.append([])
        x_sub = node_coords[ii_sub][0]
        y_sub = node_coords[ii_sub][1]
        for jj in range(H):
            nnrc_sub = nnrc[jj]
            nmrc_sup = np.prod(nmrc[jj:-1])
            x_sup = np.floor(x_sub/nnrc_sub)
            y_sup = np.floor(y_sub/nnrc_sub)
            ii_sup = nmrc_sup*x_sup+y_sup
            map_to_upper[ii_sub].append(ii_sup.astype(int))
    spatial_information['map_to_upper'] = map_to_upper      

    # index mapping between corner and degree
    index_remapping__from_corner_to_degree = np.zeros([Nh])
    for ii in range(Nh):
        index_remapping__from_corner_to_degree[ii] = np.sqrt(Nh)*node_coords[ii][0]+node_coords[ii][1] # num_row_col__nodes*node_coords[ii][0]+node_coords[ii][1]        
    index_remapping__from_degree_to_corner = np.argsort(index_remapping__from_corner_to_degree)       
    spatial_information['index_remapping__from_corner_to_degree'] = index_remapping__from_corner_to_degree
    spatial_information['index_remapping__from_degree_to_corner'] = index_remapping__from_degree_to_corner
    
    return spatial_information

def neuron_level_rentian_scaling(h,s_i,o_d_d,r_e):
        
    rent = dict()
    rent['exponent'] = r_e
    
    _tn = 0
    for ii in range(h['H__num_levels_hier']-1):
        _tn += h['n_h__num_nodes_per_module_vs_hierarchy'][ii]**( r_e-1 )                
    rent['norm_factor'] = _tn**(-1)
    # print(rent['norm_factor'])   
    
    # calculate out-degree at each level of hierarchy in accordance with rent's rule
    o_d_d['node_degrees_vs_h'] = np.zeros([h['N_h__total_num_nodes'],h['H__num_levels_hier']-1])    
    for ii in range(h['N_h__total_num_nodes']):
        for jj in range(h['H__num_levels_hier']-1):
            o_d_d['node_degrees_vs_h'][ii,jj] = np.round( o_d_d['node_degrees'][ii]*rent['norm_factor'] * ( h['n_h__num_nodes_per_module_vs_hierarchy'][jj]**( r_e-1 ) ) ).astype(int)

    # assign connections at each level of hierarchy with specified spatial profile
    A = np.zeros([h['N_h__total_num_nodes'],h['N_h__total_num_nodes']])    
    
    print('constructing adjacency matrix ...')
    for ii in range(h['N_h__total_num_nodes']):
        # print('ii = {} of {}'.format(ii+1,h['total_num_nodes']))
        
        pp = s_i['index_remapping__from_corner_to_degree'][ii].astype(int)
        # print('ii = {} of {}; pp = {}'.format(ii,h['total_num_nodes'],pp))
                    
        for jj in range(h['H__num_levels_hier']-1):
                
            module_index = s_i['map_to_upper'][ii][jj+1]
            # print('jj = {}; module_index = {}'.format(jj,module_index))
            candidate_indices = s_i['intra_modular_indices'][jj+1][module_index]
            
            # samples = np.random.randint(0,len(np.asarray(candidate_indices)),o_d_d['node_degrees_vs_h'][ii,jj].astype(int))
            # print('samples = {}'.format(samples))
            # print('candidate_indices = {}'.format(candidate_indices))
            
            num_nodes_placed = 0
            while num_nodes_placed < o_d_d['node_degrees_vs_h'][ii,jj]:
                
                sub_module_index = s_i['map_to_upper'][ii][jj]
                non_candidate_indices = s_i['intra_modular_indices'][jj][sub_module_index]
                # print('non_candidate_indices = {}'.format(non_candidate_indices))
                
                sample = np.random.randint(0,len(np.asarray(candidate_indices)),1)[0]
                
                if candidate_indices[sample].astype(int) not in np.asarray(non_candidate_indices):
                    A[pp,candidate_indices[sample].astype(int)] = 1
                    num_nodes_placed += 1
                # else:
                #     print('a node was rejected')

    return A, rent, o_d_d

def graph_analysis(A,hierarchy,s_i,rent,analyze_path_length = False,analyze_small_world = False,modularity_matrix = False):
    
    print('performing graph analysis ...')
    
    graph_data = dict()
    
    print('  initializing networkx object ...')
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    tot_num_nodes = G.number_of_nodes()
            
    if analyze_path_length == True:
        print('  calculating average shortest path length ...')
        st = time.time()
        spl = nx.average_shortest_path_length(G) # networkx calculation of the average shortest path length. # nx.average_shortest_path_length(G[, weight]) 
        print('    that took {:7.2}s'.format(time.time()-st))    
        graph_data['average_shortest_path_length'] = spl
    
    print('  analyzing node degrees ...')
    in_degree = G.in_degree()
    out_degree = G.out_degree()
    
    in_degree_vec = np.zeros([tot_num_nodes])
    out_degree_vec = np.zeros([tot_num_nodes])
    
    tot_in_degree = 0
    tot_out_degree = 0
    
    for ii in range(tot_num_nodes):
        
        in_degree_vec[ii] = in_degree[ii]
        tot_in_degree += in_degree_vec[ii]
        
        out_degree_vec[ii] = out_degree[ii]
        tot_out_degree += out_degree[ii]
        
    if tot_in_degree != tot_out_degree:
        raise ValueError('[_functions_network/graph_analysis] The total in-degree of the graph does not equal the total out-degree')
    avg_in_degree = tot_in_degree/tot_num_nodes        
    avg_out_degree = tot_out_degree/tot_num_nodes
    
    standard_deviation_in_degree = np.sqrt( np.sum( (in_degree_vec - avg_in_degree)**2 )/tot_num_nodes )
    standard_deviation_out_degree = np.sqrt( np.sum( (out_degree_vec - avg_out_degree)**2 )/tot_num_nodes )
    
    print('  generating random graph for comparison ...')
    G__rand = nx.fast_gnp_random_graph(tot_num_nodes, tot_out_degree/(tot_num_nodes*(tot_num_nodes-1)), seed = None, directed = True) 
    
    if analyze_path_length == True:
        print('  calculating random graph shortest path length ...')
        st = time.time()
        spl__rand = nx.average_shortest_path_length(G__rand)
        print('    that took {:7.2}s'.format(time.time()-st))
        graph_data['average_shortest_path_length__random_graph'] = spl__rand
    
    if analyze_small_world == True:
        print('  analyzing clustering ...')
        st = time.time()
        # clustering =  nx.clustering(G)
        average_clustering_coefficient =  nx.average_clustering(G)
        average_clustering_coefficient__rand =  nx.average_clustering(G__rand)
        print('    that took {:7.2}s for real network and random combined'.format(time.time()-st))
        graph_data['average_clustering_coefficient'] = average_clustering_coefficient
        graph_data['average_clustering_coefficient__random_graph'] = average_clustering_coefficient__rand
        
        print('  calculating small-world index ...')
        swi = (average_clustering_coefficient*spl__rand)/(average_clustering_coefficient__rand*spl)
        print('    swi = {}'.format(swi))
        graph_data['small_world_index'] = swi
    
    print('  performing rentian analysis ...')
    st = time.time()
    nlh = hierarchy['H__num_levels_hier']
    e_h_hp1 = np.zeros([nlh-1])
    # print('np.shape(e_h_hp1) = {}'.format(np.shape(e_h_hp1)))
    connections_list = []
    
    for ii in range(tot_num_nodes):
        connections = np.where( A[ii,:] )[0]
        connections_list.append(connections)
        
        for jj in range(len(connections)):
            
            h = 0
            loop_breaker = 0
            while loop_breaker == 0 and h < nlh:
                # print('connections[jj] = {}'.format(connections[jj]))
                # print('np.asarray(s_i[''intra_modular_indices''][h+1][s_i[''map_to_upper''][ii][h]]) = {}'.format(np.asarray(s_i['intra_modular_indices'][h+1][s_i['map_to_upper'][ii][h]])))
                if connections[jj] in np.asarray(s_i['intra_modular_indices'][h+1][s_i['map_to_upper'][ii][h+1]]):
                    e_h_hp1[h] += 1
                    loop_breaker = 1
                else:
                    h += 1
    
    # print('np.shape(e_h_hp1) = {}'.format(np.shape(e_h_hp1)))
    
    #fit to power law
    
    num_nodes_per_module__dense = np.linspace(hierarchy['n_h__num_nodes_per_module_vs_hierarchy'][0],hierarchy['n_h__num_nodes_per_module_vs_hierarchy'][-2],100)
    
    e_fit = np.polyfit(np.log10(hierarchy['n_h__num_nodes_per_module_vs_hierarchy'][0:-1]),np.log10(e_h_hp1),1)    
    rentian_prefactor = 10**(e_fit[1])
    rentian_exponent = e_fit[0]
    e_h_hp1__dense = rentian_prefactor*num_nodes_per_module__dense**rentian_exponent
    print('    that took {:7.2}s for rentian analysis'.format(time.time()-st))
    print('    rent exponent = {:4.2f}, targeting {:4.2f}'.format(rentian_exponent,rent['exponent']))
        
    if modularity_matrix:
        print('  calculating modularity matrix ...')
        B = nx.directed_modularity_matrix(G)
        
        graph_data['B'] = B
        g = ig.Graph.Adjacency(A.tolist())
        communities =  ig.GraphBase.community_infomap(g, edge_weights=None, vertex_weights=None,trials=10)
        graph_data['communities'] = communities

    
    # this throws an error
    # from networkx.algorithms.community import greedy_modularity_communities
    # communities = list(greedy_modularity_communities(G))
    
    # this doesnt works for digraphs
    # partition  = community_louvain.best_partition(G)
    # # draw the graph
    # pos = nx.spring_layout(G)
    # # color the nodes according to their partition
    # cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
    # nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40, cmap=cmap, node_color=list(partition.values()))
    # nx.draw_networkx_edges(G, pos, alpha=0.5)
    # plt.show()

    
    graph_data['G'] = G
    graph_data['in_degree'] = in_degree_vec
    graph_data['out_degree'] = out_degree_vec
    graph_data['tot_in_degree'] = tot_in_degree
    graph_data['tot_out_degree'] = tot_out_degree
    graph_data['avg_in_degree'] = avg_in_degree
    graph_data['avg_out_degree'] = avg_out_degree
    graph_data['standard_deviation_in_degree'] = standard_deviation_in_degree
    graph_data['standard_deviation_out_degree'] = standard_deviation_out_degree
    # graph_data['clustering'] = clustering
    
    graph_data['connections'] = connections_list
    graph_data['e_h_hp1'] = e_h_hp1
    graph_data['num_nodes_per_module__dense'] = num_nodes_per_module__dense
    graph_data['e_h_hp1__dense'] = e_h_hp1__dense
    graph_data['rentian_prefactor'] = rentian_prefactor
    graph_data['rentian_exponent'] = rentian_exponent
    
    
    
    return graph_data
