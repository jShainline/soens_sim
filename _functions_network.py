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
    
    M_h__num_modules_vs_hierarchy = np.zeros([num_levels_hier])
    n_h__num_nodes_vs_hierarchy = np.zeros([num_levels_hier])
    # n_h__num_nodes_vs_hierarchy[0] = 1
    
    num_nodes_0 = num_row_col_0**2
    a = delta_num_row_col
    b = num_row_col_0+delta_num_row_col
    
    h_vec = np.arange(1,num_levels_hier+1,1)
    # print('h_vec = {}'.format(h_vec))
    for h in h_vec:
        # M_h__num_modules_vs_hierarchy[h-1] = np.round( (11 - 2*h)**2 )
        M_h__num_modules_vs_hierarchy[h-1] = np.round( ( b-a*h )**2 )
        # print('np.prod( M_h__num_modules_vs_hierarchy[1:h] ) = {}'.format(np.prod( M_h__num_modules_vs_hierarchy[1:h] )))
        n_h__num_nodes_vs_hierarchy[h-1] = num_nodes_0*np.prod( M_h__num_modules_vs_hierarchy[1:h] )
    
    num_nodes_per_module = n_h__num_nodes_vs_hierarchy/M_h__num_modules_vs_hierarchy    
    total_nodes = n_h__num_nodes_vs_hierarchy[-1]
    inter_modular_nodes = n_h__num_nodes_vs_hierarchy*(1-1/M_h__num_modules_vs_hierarchy)
        
    nnrc = []
    nmrc = []
    for ii in range(num_levels_hier):
        nnrc.append( np.sqrt(num_nodes_per_module[ii]) ) # num nodes in each row/col within a module at the lower level of hierarchy
        nmrc.append( np.sqrt(M_h__num_modules_vs_hierarchy[ii]) ) # num modules in each row/col at this level of hierarchy
        
    # num_sub_modules_list = np.zeros([num_levels_hier-1])
    # for ii in range(num_levels_hier-1):
    #     num_sub_modules_list[ii] = M_h__num_modules_vs_hierarchy[ii]/M_h__num_modules_vs_hierarchy[ii+1]
           
    hierarchy = dict()
    
    hierarchy['num_nodes_0'] = num_nodes_0
    hierarchy['H__num_levels_hier'] = num_levels_hier
    hierarchy['h_vec'] = h_vec
    
    hierarchy['M_h__num_modules_vs_hierarchy'] = M_h__num_modules_vs_hierarchy.astype(int) # number of modules at each level of hierarchy
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
    
    # assign node coordinates within module
    # num_row_col__nodes__level_1 = np.sqrt(hierarchy['num_modules_list'][0]).astype(int)
    # num_nodes_per_module__level_1 = hierarchy['num_modules_list'][0]
    
    Nh = hier['N_h__total_num_nodes']
    H = hier['H__num_levels_hier']
    Mh = hier['M_h__num_modules_vs_hierarchy']
    nh = hier['n_h__num_nodes_vs_hierarchy']
    nmrc = hier['num_modules_row_col']
    nnrc = hier['num_nodes_row_col']
    
    # first establish modules and available sites for nodes within modules, all corner referenced to bottom left
    site_coords = []
    for ii in range(Nh):
        x = np.floor(ii/nnrc[-1])
        y = ii - nnrc[-1]*x
        site_coords.append([x,y])
        
        # x_sup = np.floor(jj/nmrc_sup)
        # y_sup = jj-nmrc_sup*x_sup
        # xind_1 = x_sup*nnrc_sub
        # xind_2 = xind_1+nnrc_sub
        # yind_1 = y_sup*nnrc_sub
        # yind_2 = yind_1+nnrc_sub
    
    spatial_information['site_coords'] = site_coords
    
    mod_center_coords = []
    mod_corner_coords = []
    mod_index = []
    for h in range(H):
        mod_center_coords.append([])
        mod_corner_coords.append([])
        for mh in range(np.prod(Mh[h:H-1])):
            x_mod = np.floor(mh/np.prod(nmrc[h:H-1]))
            y_mod = mh-np.prod(nmrc[h:H-1])*x_mod
            delta = np.round( (nnrc[h]-0.1)/2 )
            mod_center_coords[h].append([x_mod*nnrc[h]+delta,y_mod*nnrc[h]+delta])
            mod_corner_coords[h].append([x_mod*nnrc[h],y_mod*nnrc[h]])
            # mod_center_coords[h].append([x_mod*nnrc[h]+delta,y_mod*nnrc[h]+delta])
            # mod_corner_coords[h].append([x_mod,y_mod])
            
    spatial_information['mod_center_coords'] = mod_center_coords 
    spatial_information['mod_corner_coords'] = mod_corner_coords       
        
    
    # mod_center_coords = []
    # for h in range(Nh):
    #     mod_center_coords.append([])
    #     for mh in range(Mh[H-h-1]):
            
    #         # central_node_index = np.round( (np.sqrt(hier['num_modules_list'][0])-1)/2 + 1 )
    #         mod_center_coords.append([])
        
        
    # num_modules_level_1 = np.prod(hierarchy['num_modules_list'][1:-1])
    # num_row_col__modules = np.sqrt(num_modules_level_1).astype(int)
    
    # module_index__start_corner = np.zeros([num_modules_level_1])
    # module_coords__start_corner = []
    # module_coords__start_center = []
    # for ii in range(num_modules_level_1):
    #     module_index__start_corner[ii] = ii 
    #     nn = np.floor(ii/num_row_col__modules)
    #     module_coords__start_corner.append([ii-num_row_col__modules*nn,nn])
    
    # central_module_index = np.round( (np.sqrt(num_modules_level_1)-1)/2 +1 )
    # central_module_coords = [central_module_index-1,central_module_index-1]
    # module_coords_list = copy.deepcopy(module_coords__start_corner)
    # for ii in range(num_modules_level_1):
        
    #     distance_list = np.zeros([len(module_coords_list)])
    #     for jj in range(len(module_coords_list)):
    #         distance_list[jj] = ( (module_coords_list[jj][0]-central_module_coords[0])**2 + (module_coords_list[jj][1]-central_module_coords[1])**2 )**(1/2) # euclidean distance
    
    #     ind = np.argmin( distance_list )
    #     module_coords__start_center.append( module_coords_list[ind] )
    #     module_coords_list = np.delete(module_coords_list,ind,0)    
    
    # num_row_col__nodes = hier['num_row_col']
    # map_to_upper = []
    # for ii_sub in range(Nh):
        
    #     map_to_upper.append([])
        
    #     x_sub = node_coords[ii_sub][0]
    #     y_sub = node_coords[ii_sub][1]
        
    #     for jj in range(H):
    #         nnrc_sub = nnrc[jj]
    #         nmrc_sup = np.prod(nmrc[jj:-1])
    #         x_sup = np.floor(x_sub/nnrc_sub)
    #         y_sup = np.floor(y_sub/nnrc_sub)
    #         ii_sup = nmrc_sup*x_sup+y_sup
    #         map_to_upper[ii_sub].append(ii_sup.astype(int))

        

    # for ii in range(nlh):
    #     central_node_index = np.round( (np.sqrt(hier['num_modules_list'][0])-1)/2 + 1 )
        
    # coords_list = []
    # for ii in range(num_row_col__nodes__level_1):
    #     for jj in range(num_row_col__nodes__level_1):
    #         coords_list.append(np.asarray([ii,jj]))
            
    
    # c_coords = [central_node_index-1,central_node_index-1]
    # intra_module_coords__template = []
    # for ii in range(num_nodes_per_module__level_1):
        
    #     distance_list = np.zeros([len(coords_list)])
    #     for jj in range(len(coords_list)):        
    #         distance_list[jj] = ( (coords_list[jj][0]-c_coords[0])**2 + (coords_list[jj][1]-c_coords[1])**2 )**(1/2) # euclidean distance
    
    #     ind = np.argmin( distance_list )
    #     intra_module_coords__template.append( coords_list[ind] )
    #     coords_list = np.delete(coords_list,ind,0) 
    

      
    # # go through all modules, assigning nodes in descending order of out degree
    # node_coords = []    
    # num_nodes_per_module = hierarchy['num_modules_list'][0]
    # tn = num_row_col__nodes__level_1
    # ta1 = module_coords__start_center
    # ta2 = intra_module_coords__template
    # # print(np.shape(ta1))
    # # print(num_modules_level_1)
    # # print(np.shape(ta2))
    # # print(num_nodes_per_module)
    # for ii in range(num_nodes_per_module):
    #     for jj in range(num_modules_level_1):
    #         node_coords.append( [ tn*ta1[jj][0]+ta2[ii][0] , tn*ta1[jj][1]+ta2[ii][1]  ] ) 
            
             

            
    # intra_modular_indices = []   
    # # num_modules_list__with_level_0 = np.insert(copy.deepcopy(hierarchy['num_modules_list']),0,tnn)
    # # print('num_modules_list__with_level_0 = {}'.format(num_modules_list__with_level_0))
    # module_coords__x = []
    # module_coords__y = []    
    # h_increment = 0.2
    # for ii in range(nlh):
    #     intra_modular_indices.append([])
    #     module_coords__x.append([])
    #     module_coords__y.append([])
        
    #     if ii == 0:
    #         for jj in range(tnn):
    #             intra_modular_indices[ii].append(jj)
    #             module_coords__x[ii].append(node_coords[jj][0])
    #             module_coords__y[ii].append(node_coords[jj][0])
    #     elif ii > 0:
                
    #         nnrc_sub = hierarchy['num_nodes_row_col'][ii]
    #         nmrc_sup = np.prod(hierarchy['num_modules_row_col'][ii:])
            
    #         # for jj in range(num_modules_list__with_level_0):
    #         for jj in range(np.prod(hierarchy['num_modules_list'][ii:])):
                
    #             x_sup = np.floor(jj/nmrc_sup)
    #             y_sup = jj-nmrc_sup*x_sup
    #             xind_1 = x_sup*nnrc_sub
    #             xind_2 = xind_1+nnrc_sub
    #             yind_1 = y_sup*nnrc_sub
    #             yind_2 = yind_1+nnrc_sub
                
    #             # print('xind_1 = {}, xind_2 = {}, yind_1 = {}, yind_2 = {}'.format(xind_1,xind_2,yind_1,yind_2))
                
    #             xvec = np.arange(xind_1,xind_2,1)
    #             yvec = np.arange(yind_1,yind_2,1)
    #             # temp_vec = np.zeros([( (xind_2-xind_1)*(yind_2-yind_1) ).astype(int)])
    #             module_coords__x[ii].append([xind_1-ii*h_increment,xind_2-1+ii*h_increment,xind_2-1+ii*h_increment,xind_1-ii*h_increment,xind_1-ii*h_increment])
    #             module_coords__y[ii].append([yind_1-ii*h_increment,yind_1-ii*h_increment,yind_2-1+ii*h_increment,yind_2-1+ii*h_increment,yind_1-ii*h_increment])
    #             temp_vec = []
    #             for pp in range(len(xvec)):
    #                 for qq in range(len(yvec)):
    #                     # temp_vec[pp*len(yvec)+qq] = nnrc_sub*xvec[pp]+yvec[qq]
    #                     temp_vec.append(num_row_col__nodes*xvec[pp]+yvec[qq])
                
    #             intra_modular_indices[ii].append(temp_vec)
            
    # inter_modular_indices = [] 
    # # print('removing nodes from inter_modular_indices ...')
    # # for ii in range(nlh-1):
    # #     inter_modular_indices.append(intra_modular_indices[ii+1]) 
    # #     for jj in range(nlh-1):
    # #         print('ii = {}, jj = {}, map_to_upper[ii][jj] = {}'.format(ii,jj,map_to_upper[ii][jj]))
    # #         index = np.where( inter_modular_indices[jj][map_to_upper[ii][jj]] == ii )
    # #         inter_modular_indices[jj][map_to_upper[ii][jj]] = np.delete(inter_modular_indices[jj][map_to_upper[ii][jj]],index)
    
                        
    # # print(np.shape(node_coords))
    # degree_xy = np.zeros([num_row_col__nodes,num_row_col__nodes])
    # for ii in range(tnn):
    #     # print('mm_ii = {}, nn_ii = {}'.format(mm_ii,nn_ii))
    #     degree_xy[node_coords[ii][0].astype(int),node_coords[ii][1].astype(int)] = out_degree_distribution['node_degrees'][ii]
                        
    # # print(node_coords[0])
    # index_remapping__from_corner_to_degree = np.zeros([tnn])
    # for ii in range(tnn):
    #     index_remapping__from_corner_to_degree[ii] = num_row_col__nodes*node_coords[ii][0]+node_coords[ii][1]
        
    # index_remapping__from_degree_to_corner = np.argsort(index_remapping__from_corner_to_degree)
    
    # # generate distance matrix

    # # print('generating distance matrix ... ')
    
    # # distance_mat = np.zeros([tnn,tnn])
    # # node_coords__corner = []
    # # for ii in range(tnn):
    # #     # print('ii = {} of {} (num_nodes)'.format(ii+1,total_num_nodes))
    # #     mm = np.floor(ii/hierarchy['num_nodes_row_col'][-1])
    # #     node_coords__corner.append( [mm,ii-hierarchy['num_nodes_row_col'][-1]*mm] )
        
    # #     for jj in range(tnn):
    # #         distance_mat[ii,jj] = np.abs(node_coords[ii][0]-node_coords[jj][0])+np.abs(node_coords[ii][1]-node_coords[jj][1]) # manhattan distance
    # #         # distance_mat[ii,jj] = ( (mm_ii[ii]-mm_ii[jj])**2 + (nn_ii[ii]-nn_ii[jj])**2 )**(1/2) # euclidean distance            
            
    # # distance_mat__corner = np.zeros([tnn,tnn])
    # # for ii in range(tnn):
    # #     # print('ii = {} of {} (num_nodes)'.format(ii+1,total_num_nodes))
    # #     for jj in range(tnn):
    # #         # distance_mat__corner[ii,jj] = np.abs(node_coords__corner[ii][0]-node_coords__corner[jj][0])+np.abs(node_coords__corner[ii][1]-node_coords__corner[jj][1]) # manhattan distance
    # #         distance_mat__corner[ii,jj] = distance_mat[index_remapping__from_degree_to_corner[ii].astype(int),index_remapping__from_degree_to_corner[jj].astype(int)]

    
    # spatial_information['mod_center_coords'] = mod_center_coords
    # spatial_information['mod_node_coords'] = mod_node_coords
    # # spatial_information['module_index__start_corner'] = module_index__start_corner
    # spatial_information['module_coords__x'] = module_coords__x
    # spatial_information['module_coords__y'] = module_coords__y
    # # spatial_information['module_coords__start_center'] = module_coords__start_center
    # spatial_information['node_coords'] = node_coords
    # spatial_information['degree_xy'] = degree_xy
    # # spatial_information['distance_mat'] = distance_mat
    # # spatial_information['distance_mat__corner'] = distance_mat__corner
    # spatial_information['map_to_upper'] = map_to_upper
    # spatial_information['intra_modular_indices'] = intra_modular_indices
    # spatial_information['inter_modular_indices'] = inter_modular_indices
    # spatial_information['index_remapping__from_degree_to_corner'] = index_remapping__from_degree_to_corner.astype(int)
    # spatial_information['index_remapping__from_corner_to_degree'] = index_remapping__from_corner_to_degree.astype(int)
    
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


def neuron_level_rentian_scaling(h,s_i,o_d_d,r_e):
        
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


def graph_analysis(A,hierarchy,s_i,rent):
    
    print('performing graph analysis ...')
    
    print('  initializing networkx object ...')
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    tot_num_nodes = G.number_of_nodes()
            
    print('  calculating average shortest path length ...')
    st = time.time()
    spl = nx.average_shortest_path_length(G) # networkx calculation of the average shortest path length. # nx.average_shortest_path_length(G[, weight]) 
    print('    that took {:7.2}s'.format(time.time()-st))    
    
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
    
    print('  calculating random graph shortest path length ...')
    st = time.time()
    spl__rand = nx.average_shortest_path_length(G__rand)
    print('    that took {:7.2}s'.format(time.time()-st))
    
    print('  analyzing clustering ...')
    st = time.time()
    # clustering =  nx.clustering(G)
    average_clustering_coefficient =  nx.average_clustering(G)
    average_clustering_coefficient__rand =  nx.average_clustering(G__rand)
    print('    that took {:7.2}s for real network and random combined'.format(time.time()-st))
    
    print('  calculating small-world index ...')
    swi = (average_clustering_coefficient*spl__rand)/(average_clustering_coefficient__rand*spl)
    print('    swi = {}'.format(swi))
    
    print('  performing rentian analysis ...')
    nlh = hierarchy['num_levels_hier']
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
    
    num_nodes_per_module__dense = np.linspace(hierarchy['num_nodes_per_module'][0],hierarchy['num_nodes_per_module'][-2],100)
    
    e_fit = np.polyfit(np.log10(hierarchy['num_nodes_per_module'][0:-1]),np.log10(e_h_hp1),1)    
    rentian_prefactor = 10**(e_fit[1])
    rentian_exponent = e_fit[0]
    e_h_hp1__dense = rentian_prefactor*num_nodes_per_module__dense**rentian_exponent
    print('    rent exponent = {:4.2f}, targeting {:4.2f}'.format(rentian_exponent,rent['exponent']))
        
    print('  calculating modularity matrix ...')
    B = nx.directed_modularity_matrix(G)
    
    g = ig.Graph.Adjacency(A.tolist())
    communities =  ig.GraphBase.community_infomap(g, edge_weights=None, vertex_weights=None,trials=10)

    
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

    
    
    graph_data = dict()
    graph_data['G'] = G
    graph_data['B'] = B
    graph_data['average_shortest_path_length'] = spl
    graph_data['average_shortest_path_length__random_graph'] = spl__rand
    graph_data['in_degree'] = in_degree_vec
    graph_data['out_degree'] = out_degree_vec
    graph_data['tot_in_degree'] = tot_in_degree
    graph_data['tot_out_degree'] = tot_out_degree
    graph_data['avg_in_degree'] = avg_in_degree
    graph_data['avg_out_degree'] = avg_out_degree
    graph_data['standard_deviation_in_degree'] = standard_deviation_in_degree
    graph_data['standard_deviation_out_degree'] = standard_deviation_out_degree
    # graph_data['clustering'] = clustering
    graph_data['average_clustering_coefficient'] = average_clustering_coefficient
    graph_data['average_clustering_coefficient__random_graph'] = average_clustering_coefficient__rand
    graph_data['small_world_index'] = swi
    graph_data['communities'] = communities
    graph_data['connections'] = connections_list
    graph_data['e_h_hp1'] = e_h_hp1
    graph_data['num_nodes_per_module__dense'] = num_nodes_per_module__dense
    graph_data['e_h_hp1__dense'] = e_h_hp1__dense
    graph_data['rentian_prefactor'] = rentian_prefactor
    graph_data['rentian_exponent'] = rentian_exponent
    
    
    
    return graph_data


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