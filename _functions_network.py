import numpy as np
import copy
# from matplotlib import pyplot as plt
import powerlaw

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
    num_nodes_list = np.zeros([num_levels_hier])
    num_nodes_list[0] = 1
    
    h_vec = np.arange(1,num_levels_hier+1,1)
    for h in h_vec:
        num_modules_list[h-1] = np.round( num_nodes_0 * h**(-gamma) ) # num_nodes_0 * h**(-gamma) # 
        num_nodes_list[h-1] = ( num_nodes_0**h )*( np.prod(np.arange(1,h+1,1)**(-gamma)) ) # num_nodes_0*( np.prod(num_nodes_0*np.arange(2,h+1,1)**(-gamma)) ) # 
    
    num_nodes_per_module = num_nodes_list/num_modules_list    
    total_nodes = num_nodes_list[-1]
    inter_modular_nodes = num_nodes_list*(1-1/num_modules_list)

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
    hierarchy['num_nodes_list'] = num_nodes_list # number of nodes at each level of hierarchy
    hierarchy['num_nodes_per_module'] = num_nodes_per_module # number of nodes per module at each level of hierarchy
    hierarchy['inter_modular_nodes'] = inter_modular_nodes # number of inter-modular nodes at each level of hierarchy
    hierarchy['total_num_nodes'] = total_nodes # total number of nodes in the network 
    hierarchy['num_nodes_row_col'] = nnrc
    hierarchy['num_modules_row_col'] = nmrc        

    return hierarchy


def populate_hierarchy__geometrical(num_row_col_0 = 9, delta_num_row_col = 2, num_levels_hier = 4, plot = True): # power-law hierarchy ( number of modules at level h of hierarchy: M_h = n_0 * h**(-gamma) )
    
    num_modules_list = np.zeros([num_levels_hier])
    num_nodes_list = np.zeros([num_levels_hier])
    # num_nodes_list[0] = 1
    
    num_nodes_0 = num_row_col_0**2
    a = delta_num_row_col
    b = num_row_col_0+delta_num_row_col
    
    h_vec = np.arange(1,num_levels_hier+1,1)
    # print('h_vec = {}'.format(h_vec))
    for h in h_vec:
        # num_modules_list[h-1] = np.round( (11 - 2*h)**2 )
        num_modules_list[h-1] = np.round( ( b-a*h )**2 )
        # print('np.prod( num_modules_list[1:h] ) = {}'.format(np.prod( num_modules_list[1:h] )))
        num_nodes_list[h-1] = num_nodes_0*np.prod( num_modules_list[1:h] )
    
    num_nodes_per_module = num_nodes_list/num_modules_list    
    total_nodes = num_nodes_list[-1]
    inter_modular_nodes = num_nodes_list*(1-1/num_modules_list)
        
    nnrc = []
    nmrc = []
    for ii in range(num_levels_hier):
        nnrc.append( np.sqrt(num_nodes_per_module[ii]) ) # num nodes in each row/col within a module at the lower level of hierarchy
        nmrc.append( np.sqrt(num_modules_list[ii]) ) # num modules in each row/col at this level of hierarchy
        
    # num_sub_modules_list = np.zeros([num_levels_hier-1])
    # for ii in range(num_levels_hier-1):
    #     num_sub_modules_list[ii] = num_modules_list[ii]/num_modules_list[ii+1]
           
    hierarchy = dict()
    
    hierarchy['num_nodes_0'] = num_nodes_0
    hierarchy['num_levels_hier'] = num_levels_hier
    hierarchy['h_vec'] = h_vec
    
    hierarchy['num_modules_list'] = num_modules_list.astype(int) # number of modules at each level of hierarchy
    # hierarchy['num_sub_modules_list'] = num_sub_modules_list.astype(int) # number of sub-modules contained in a module at each level of hierarchy above the lowest
    # hierarchy['map_to_upper'] = map_to_upper.astype(int) # index of containing module next level up in hierachy for each level of hierarchy except the highest
    hierarchy['num_nodes_list'] = num_nodes_list.astype(int) # number of nodes at each level of hierarchy
    hierarchy['num_nodes_per_module'] = num_nodes_per_module.astype(int) # number of nodes per module at each level of hierarchy
    hierarchy['inter_modular_nodes'] = inter_modular_nodes.astype(int) # number of inter-modular nodes at each level of hierarchy
    hierarchy['total_num_nodes'] = total_nodes.astype(int) # total number of nodes in the network         
    hierarchy['num_row_col'] = np.sqrt(total_nodes).astype(int) # number of rows and columns in the square network    
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


def generate_spatial_structure(hierarchy,out_degree_distribution):

    # assign node coordinates within module
    num_row_col__nodes__level_1 = np.sqrt(hierarchy['num_modules_list'][0]).astype(int)
    num_nodes_per_module__level_1 = hierarchy['num_modules_list'][0]
        
    coords_list = []
    for ii in range(num_row_col__nodes__level_1):
        for jj in range(num_row_col__nodes__level_1):
            coords_list.append(np.asarray([ii,jj]))
            
    central_node_index = np.round( (np.sqrt(hierarchy['num_modules_list'][0])-1)/2 +1 )
    c_coords = [central_node_index-1,central_node_index-1]
    intra_module_coords__template = []
    for ii in range(num_nodes_per_module__level_1):
        
        distance_list = np.zeros([len(coords_list)])
        for jj in range(len(coords_list)):        
            distance_list[jj] = ( (coords_list[jj][0]-c_coords[0])**2 + (coords_list[jj][1]-c_coords[1])**2 )**(1/2) # euclidean distance
    
        ind = np.argmin( distance_list )
        intra_module_coords__template.append( coords_list[ind] )
        coords_list = np.delete(coords_list,ind,0) 
    
    # assign module coordinates relative to bottom left
    num_modules_level_1 = np.prod(hierarchy['num_modules_list'][1:-1])
    num_row_col__modules = np.sqrt(num_modules_level_1).astype(int)
    
    module_index__start_corner = np.zeros([num_modules_level_1])
    module_coords__start_corner = []
    module_coords__start_center = []
    for ii in range(num_modules_level_1):
        module_index__start_corner[ii] = ii 
        nn = np.floor(ii/num_row_col__modules)
        module_coords__start_corner.append([ii-num_row_col__modules*nn,nn])
    
    central_module_index = np.round( (np.sqrt(num_modules_level_1)-1)/2 +1 )
    central_module_coords = [central_module_index-1,central_module_index-1]
    module_coords_list = copy.deepcopy(module_coords__start_corner)
    for ii in range(num_modules_level_1):
        
        distance_list = np.zeros([len(module_coords_list)])
        for jj in range(len(module_coords_list)):
            distance_list[jj] = ( (module_coords_list[jj][0]-central_module_coords[0])**2 + (module_coords_list[jj][1]-central_module_coords[1])**2 )**(1/2) # euclidean distance
    
        ind = np.argmin( distance_list )
        module_coords__start_center.append( module_coords_list[ind] )
        module_coords_list = np.delete(module_coords_list,ind,0)
      
    # go through all modules, assigning nodes in descending order of out degree
    node_coords = []    
    num_nodes_per_module = hierarchy['num_modules_list'][0]
    tn = num_row_col__nodes__level_1
    ta1 = module_coords__start_center
    ta2 = intra_module_coords__template
    # print(np.shape(ta1))
    # print(num_modules_level_1)
    # print(np.shape(ta2))
    # print(num_nodes_per_module)
    for ii in range(num_nodes_per_module):
        for jj in range(num_modules_level_1):
            node_coords.append( [ tn*ta1[jj][0]+ta2[ii][0] , tn*ta1[jj][1]+ta2[ii][1]  ] ) 
            
             
    tnn = hierarchy['total_num_nodes']
    nlh = hierarchy['num_levels_hier']
    num_row_col__nodes = hierarchy['num_row_col']
    map_to_upper = []
    for ii_sub in range(tnn):
        map_to_upper.append([])
        
        x_sub = node_coords[ii_sub][0]
        y_sub = node_coords[ii_sub][1]
        for jj in range(nlh-1):
            nnrc_sub = hierarchy['num_nodes_row_col'][jj+1]
            nmrc_sup = np.prod(hierarchy['num_modules_row_col'][jj+1:-1])
            x_sup = np.floor(x_sub/nnrc_sub)
            y_sup = np.floor(y_sub/nnrc_sub)
            ii_sup = nmrc_sup*x_sup+y_sup
            map_to_upper[ii_sub].append(ii_sup)
            
    intra_modular_indices = []   
    # num_modules_list__with_level_0 = np.insert(copy.deepcopy(hierarchy['num_modules_list']),0,tnn)
    # print('num_modules_list__with_level_0 = {}'.format(num_modules_list__with_level_0))
    for ii in range(nlh-1):
        intra_modular_indices.append([])
        
        nnrc_sub = hierarchy['num_nodes_row_col'][ii+1]
        nmrc_sup = np.prod(hierarchy['num_modules_row_col'][ii+1:-1])
        
        # for jj in range(num_modules_list__with_level_0):
        for jj in range(np.prod(hierarchy['num_modules_list'][ii+1:-1])):
            
            x_sup = np.floor(jj/nmrc_sup)
            y_sup = jj-nmrc_sup*x_sup
            xind_1 = x_sup*nnrc_sub
            xind_2 = xind_1+nnrc_sub
            yind_1 = y_sup*nnrc_sub
            yind_2 = yind_1+nnrc_sub
            
            # print('xind_1 = {}, xind_2 = {}, yind_1 = {}, yind_2 = {}'.format(xind_1,xind_2,yind_1,yind_2))
            
            xvec = np.arange(xind_1,xind_2,1)
            yvec = np.arange(yind_1,yind_2,1)
            # temp_vec = np.zeros([( (xind_2-xind_1)*(yind_2-yind_1) ).astype(int)])
            temp_vec = []
            for pp in range(len(xvec)):
                for qq in range(len(yvec)):
                    # temp_vec[pp*len(yvec)+qq] = nnrc_sub*xvec[pp]+yvec[qq]
                    temp_vec.append(num_row_col__nodes*xvec[pp]+yvec[qq])
            
            intra_modular_indices[ii].append(temp_vec)
                        
    # print(np.shape(node_coords))
    degree_xy = np.zeros([num_row_col__nodes,num_row_col__nodes])
    for ii in range(tnn):
        # print('mm_ii = {}, nn_ii = {}'.format(mm_ii,nn_ii))
        degree_xy[node_coords[ii][0].astype(int),node_coords[ii][1].astype(int)] = out_degree_distribution['node_degrees'][ii]
                        
    # print(node_coords[0])
    index_remapping__from_corner_to_degree = np.zeros([tnn])
    for ii in range(tnn):
        index_remapping__from_corner_to_degree[ii] = num_row_col__nodes*node_coords[ii][0]+node_coords[ii][1]
        
    index_remapping__from_degree_to_corner = np.argsort(index_remapping__from_corner_to_degree)
    
    # generate distance matrix

    # print('generating distance matrix ... ')
    
    # distance_mat = np.zeros([tnn,tnn])
    # node_coords__corner = []
    # for ii in range(tnn):
    #     # print('ii = {} of {} (num_nodes)'.format(ii+1,total_num_nodes))
    #     mm = np.floor(ii/hierarchy['num_nodes_row_col'][-1])
    #     node_coords__corner.append( [mm,ii-hierarchy['num_nodes_row_col'][-1]*mm] )
        
    #     for jj in range(tnn):
    #         distance_mat[ii,jj] = np.abs(node_coords[ii][0]-node_coords[jj][0])+np.abs(node_coords[ii][1]-node_coords[jj][1]) # manhattan distance
    #         # distance_mat[ii,jj] = ( (mm_ii[ii]-mm_ii[jj])**2 + (nn_ii[ii]-nn_ii[jj])**2 )**(1/2) # euclidean distance            
            
    # distance_mat__corner = np.zeros([tnn,tnn])
    # for ii in range(tnn):
    #     # print('ii = {} of {} (num_nodes)'.format(ii+1,total_num_nodes))
    #     for jj in range(tnn):
    #         # distance_mat__corner[ii,jj] = np.abs(node_coords__corner[ii][0]-node_coords__corner[jj][0])+np.abs(node_coords__corner[ii][1]-node_coords__corner[jj][1]) # manhattan distance
    #         distance_mat__corner[ii,jj] = distance_mat[index_remapping__from_degree_to_corner[ii].astype(int),index_remapping__from_degree_to_corner[jj].astype(int)]

    spatial_information = dict()
    # spatial_information['module_index__start_corner'] = module_index__start_corner
    # spatial_information['module_coords__start_corner'] = module_coords__start_corner
    # spatial_information['module_coords__start_center'] = module_coords__start_center
    spatial_information['node_coords'] = node_coords
    spatial_information['degree_xy'] = degree_xy
    # spatial_information['distance_mat'] = distance_mat
    # spatial_information['distance_mat__corner'] = distance_mat__corner
    spatial_information['map_to_upper'] = map_to_upper
    spatial_information['intra_modular_indices'] = intra_modular_indices
    spatial_information['index_remapping__from_degree_to_corner'] = index_remapping__from_degree_to_corner.astype(int)
    spatial_information['index_remapping__from_corner_to_degree'] = index_remapping__from_corner_to_degree.astype(int)
    
    return spatial_information


def determine_indices(h,si):
    
    nrc = h['num_row_col']
    
   
    # for ii in range(h['num_levels_hier']-1):
        
    #     module_max_min_coords__start_corner.append([])
    #     nnrc = h['num_nodes_row_col'][ii+1]
    #     nmrc = h['num_modules_row_col'][ii+1]
        
    #     for jj in range(h['num_modules_list'][ii+1]):            
            
    #         nn_jj = np.floor(jj/nmrc)
    #         module_max_min_coords__start_corner[ii].append( [ [nn_jj*nnrc,nn_jj*nnrc+nnrc] , [(jj-nmrc*nn_jj)*nnrc,(jj-nmrc*nn_jj)*nnrc+nnrc] ] ) # [[xmin,xmax],[ymin,ymax]]

    # make indices_arrays__intra
    indices_arrays__intra = []
    temp_vec = np.arange(0,nrc**2,1)
    map_array = np.reshape(temp_vec,[nrc,nrc])
    for ii in range(h['total_num_nodes']):
                
        indices_arrays__intra.append([])
        # indices_arrays__inter.append([])
        
        _xy = si['node_coords'][ii]
        
        for jj in range(h['num_levels_hier']):
            
            nrc = h['num_nodes_row_col'][jj]
                
            _x1 = ( np.floor(_xy[0]/nrc)*nrc ).astype(int)
            _y1 = ( np.floor(_xy[1]/nrc)*nrc ).astype(int)
            _x2 = _x1+nrc # ( np.ceil(_xy[0]/nrc)*nrc ).astype(int)
            _y2 = _y1+nrc # ( np.ceil(_xy[1]/nrc)*nrc ).astype(int)
            
            if _x2 == _x1: _x2 += 1
            if _y2 == _y1: _y2 += 1
            indices_arrays__intra[ii].append( np.concatenate(map_array[_x1:_x2,_y1:_y2]) )  

    indices_arrays__inter = []
    for ii in range(h['total_num_nodes']):
        print('ii = {} of {} (num_nodes)'.format(ii+1,h['total_num_nodes']))
                
        indices_arrays__inter.append([])
        
        for jj in range(h['num_levels_hier']-1):
                        
            indices_arrays__inter[ii].append( np.delete(indices_arrays__intra[ii][jj+1], np.where( np.in1d(indices_arrays__intra[ii][jj+1],indices_arrays__intra[ii][jj]) )[0] ) )

    indices_arrays = dict()
    indices_arrays['intra'] = indices_arrays__intra
    indices_arrays['inter'] = indices_arrays__inter
    indices_arrays['map_array'] = map_array
    # indices_arrays['module_max_min_coords__start_corner'] = module_max_min_coords__start_corner # at each level of hierarcy, for each module, this array gives x-y coordinate ranges in the form [level_of_hierarchy][module_number][[xmin,xmax],[ymin,ymax]]
    
    return indices_arrays, h


def neuron_level_rentian_scaling__with_spatial_dependence(h,s_i,o_d_d,r_e):
        
    rent = dict()
    rent['exponent'] = r_e
    
    _tn = 0
    for ii in range(h['num_levels_hier']-1):
        _tn += h['num_nodes_per_module'][ii]**( r_e-1 )                
    rent['norm_factor'] = _tn**(-1)
    # print(rent['norm_factor'])   
    
    # calculate out-degree at each level of hierarchy in accordance with rent's rule
    o_d_d['node_degrees_vs_h'] = np.zeros([h['total_num_nodes'],h['num_levels_hier']-1])    
    for ii in range(h['total_num_nodes']):
        for jj in range(h['num_levels_hier']-1):
            o_d_d['node_degrees_vs_h'][ii,jj] = np.round( o_d_d['node_degrees'][ii]*rent['norm_factor'] * ( h['num_nodes_per_module'][jj]**( r_e-1 ) ) ).astype(int)

    # assign connections at each level of hierarchy with specified spatial profile
    A = np.zeros([h['total_num_nodes'],h['total_num_nodes']])    
    
    print('constructing adjacency matrix ...')
    for ii in range(h['total_num_nodes']):
        # print('ii = {} of {}'.format(ii+1,h['total_num_nodes']))
        
        pp = s_i['index_remapping__from_corner_to_degree'][ii].astype(int)
        # print('ii = {} of {}; pp = {}'.format(ii,h['total_num_nodes'],pp))
        
            
        for jj in range(h['num_levels_hier']-1):
                
            module_index = s_i['map_to_upper'][ii][jj]
            # print('jj = {}; module_index = {}'.format(jj,module_index))
            candidate_indices = s_i['intra_modular_indices'][jj][module_index.astype(int)]
            
            samples = np.random.randint(0,len(candidate_indices),o_d_d['node_degrees_vs_h'][ii,jj].astype(int))
            # print('samples = {}'.format(samples))
            # print('candidate_indices[samples] = {}'.format(candidate_indices[samples]))
            
            for kk in range(len(samples)):
                A[pp,candidate_indices[samples[kk]].astype(int)] = 1

    return A, rent, o_d_d


#%% scraps



# proximity_factor = 1
# A = np.zeros([num_nodes,num_nodes])
# for ii in range(num_nodes):
#     k_out_ii = node_degrees[ii].astype(int)
#     r_out_ii__vec = np.random.exponential(decay_length,k_out_ii) # exponential spatial decay
#     # print('ii = {} of {}, k_out_ii = {}, len(r_out_ii__vec) = {}'.format(ii+1,num_nodes,k_out_ii,len(r_out_ii__vec)))
#     for r_out_ii in r_out_ii__vec:     
#         tracker = 0
#         candidate_nodes = np.where( np.abs( R_mat[ii,:] - r_out_ii ) <= proximity_factor  )[0]
#         # print('here0')
#         if len(candidate_nodes) > 0:
#             while tracker == 0:
                
#                 # print('len(candidate_nodes) = {}'.format(len(candidate_nodes)))
#                 rand_ind = np.random.randint(0,len(candidate_nodes),1)
#                 # print('candidate_nodes = {}, rand_ind[0] = {}, candidate_nodes[rand_ind[0]] = {}'.format(candidate_nodes,rand_ind[0],candidate_nodes[rand_ind[0]]))
    
#                 if A[ii,candidate_nodes[rand_ind[0]]] == 0:
#                     A[ii,candidate_nodes[rand_ind[0]]] = 1
#                     tracker = 1
#                     # print('here1')
#                 elif A[ii,candidate_nodes[rand_ind[0]]] == 1:
#                     candidate_nodes = np.delete(candidate_nodes,rand_ind[0])
#                     # print('here2')
#                     if len(candidate_nodes) == 0:
#                         tracker = 1

# plot_A(A)