#%%
import numpy as np
from matplotlib import pyplot as plt

from _util import color_dictionary, physical_constants
colors = color_dictionary()
p = physical_constants()

# plt.close('all')

#%%

gamma = 0.5772
NVec = np.logspace(3,8,1000)
APL_vec = [2,3,4,5] # vector of average path lengths to consider
kExpectation_vec = [10,100,1000,10000]
# diameter_vec = [2 3 4 5];%vector of network diameters to consider

#%% Erdos-Renyi APL

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


# %% Erdos-Renyi diameter
# % kExpectation_mat2 = zeros(length(NVec),length(diameter_vec));
# % 
# % for ii = 1:length(diameter_vec)
# % 
# %     kExpectation_mat2(:,ii) = exp( log(NVec(:)) / diameter_vec(ii) );
# %     
# % end
# % 
# % figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
# % loglog(NVec(:),kExpectation_mat2(:,1),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
# % hold on
# % loglog(NVec(:),kExpectation_mat2(:,2),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
# % loglog(NVec(:),kExpectation_mat2(:,3),'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
# % loglog(NVec(:),kExpectation_mat2(:,4),'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
# % xlabel('N_{tot}' ,'FontSize',fontSize,'FontName','Times')
# % ylabel('<k>','FontSize',fontSize,'FontName','Times')
# % set(gca,'FontSize',fontSize,'FontName',fontName)
# % legend(sprintf('d = %g',diameter_vec(1)),sprintf('d = %g',diameter_vec(2)),sprintf('d = %g',diameter_vec(3)),sprintf('d = %g',diameter_vec(4)))
# % % title(sprintf('kMax vs number of nodes in network, gamma = %g',gamma),'FontSize',fontSize,'FontName',fontName)
# % % xlim([1e3 1e8])
# % % ylim([1e-6 1e2])
# % grid on

# %% Scale-free APL
# m = 2;
# kExpectation_mat = zeros(length(NVec),length(APL_vec));

# for ii = 1:length(APL_vec)

#     kExpectation_mat(:,ii) = exp( (log(NVec(:))+(1.5-APL_vec(ii))*log(log(NVec(:)))-1-gamma)/(APL_vec(ii)-0.5) );
    
# end

# figure('OuterPosition',[0 0 scrsz(3) scrsz(4)]);
# loglog(NVec(:),kExpectation_mat(:,1),'Color',bRGY(3,:),'LineStyle','-','LineWidth',3)
# hold on
# loglog(NVec(:),kExpectation_mat(:,2),'Color',bRGY(8,:),'LineStyle','-','LineWidth',3)
# loglog(NVec(:),kExpectation_mat(:,3),'Color',bRGY(13,:),'LineStyle','-','LineWidth',3)
# loglog(NVec(:),kExpectation_mat(:,4),'Color',bRGY(18,:),'LineStyle','-','LineWidth',3)
# xlabel('N_{tot}' ,'FontSize',fontSize,'FontName','Times')
# ylabel('<k>','FontSize',fontSize,'FontName','Times')
# set(gca,'FontSize',fontSize,'FontName',fontName)
# legend(sprintf('APL = %g',APL_vec(1)),sprintf('APL = %g',APL_vec(2)),sprintf('APL = %g',APL_vec(3)),sprintf('APL = %g',APL_vec(4)))
# title('Average node degree versus network size for a given APL, scale-free','FontSize',fontSize,'FontName',fontName)
# % xlim([1e3 1e8])
# % ylim([1e-6 1e2])
# grid on
