from pylab import *
import numpy as np

r_300 = 0.150#2.5*2.5/100/2;#
xmin_plot = 1e2
xmax_plot = 1e5
ymin_plot = 1e2
ymax_plot = 1e8

#--- varying numPlanes ---#
w = 1.5e-6

numPlanes_vec = np.array([2,6,18])
numConnections_vec = np.logspace(2,5,1000)
#numConnections_vec = np.linspace(700,80000,1000)

numNeurons_mat = np.zeros([len(numPlanes_vec),len(numConnections_vec)])
for ii in range(len(numPlanes_vec)):
    for jj in range(len(numConnections_vec)):
#        numNeurons_mat[ii,jj] = (2*np.sqrt(2)*numPlanes_vec[ii]**2*r_300**2/(w**2))/(numConnections_vec[jj]**2);
        numNeurons_mat[ii,jj] = 2*np.sqrt(2)*r_300**2*(numPlanes_vec[ii]/(w*numConnections_vec[jj]))**2;

fig, axes = plt.subplots(1,1)
axes.loglog(numConnections_vec[:],numNeurons_mat[0,:],label = '2 planes')
axes.loglog(numConnections_vec[:],numNeurons_mat[1,:],label = '6 planes')
axes.loglog(numConnections_vec[:],numNeurons_mat[2,:],label = '18 planes')
axes.loglog(numConnections_vec[:],numConnections_vec[:],label = 'numNeuron = numConn')
axes.loglog(numConnections_vec[:],10*numConnections_vec[:],label = 'numNeuron = 10 numConn')
axes.loglog(numConnections_vec[:],100*numConnections_vec[:],label = 'numNeuron = 100 numConn')
axes.legend()
axes.set_xlabel(r'connections per neuron', fontsize=20)
axes.set_ylabel(r'num neurons on 300 mm wafer', fontsize=20);
ylim((ymin_plot,ymax_plot))
xlim((xmin_plot,xmax_plot))
axes.legend(loc='best')
grid(True,which='both')
show() 


#--- varying w ---#
numPlanes= 6

w_vec = np.array([1.5e-6,3e-6,6e-6])

numNeurons_mat = np.zeros([len(w_vec),len(numConnections_vec)])
for ii in range(len(w_vec)):
    for jj in range(len(numConnections_vec)):
        numNeurons_mat[ii,jj] = 2*np.sqrt(2)*r_300**2*(numPlanes/(w_vec[ii]*numConnections_vec[jj]))**2;

fig, axes = plt.subplots(1,1)
axes.loglog(numConnections_vec[:],numNeurons_mat[0,:],label = 'w = 1.5 um')
axes.loglog(numConnections_vec[:],numNeurons_mat[1,:],label = 'w = 3 um')
axes.loglog(numConnections_vec[:],numNeurons_mat[2,:],label = 'w = 6 um')
axes.loglog(numConnections_vec[:],numConnections_vec[:],label = 'numNeuron = numConn')
axes.loglog(numConnections_vec[:],10*numConnections_vec[:],label = 'numNeuron = 10 numConn')
axes.loglog(numConnections_vec[:],100*numConnections_vec[:],label = 'numNeuron = 100 numConn')
axes.legend()
axes.set_xlabel(r'connections per neuron', fontsize=20)
axes.set_ylabel(r'num neurons on 300 mm wafer', fontsize=20);
ylim((ymin_plot,ymax_plot))
xlim((xmin_plot,xmax_plot))
axes.legend(loc='best')
grid(True,which='both')
show() 
