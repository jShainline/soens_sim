import numpy as np

from _util import physical_constants
p = physical_constants()
Phi0 = p['Phi0__pH_ns']

import WRSpice


#%%

# tempCirFile only netlist
# export your netlist to a cir file using 'deck' command in xic
# open the cir file and remove all lines except the netlist lines (include .model)
# replace all values by Params, you want to sweep, with white space around
# Component names don't need any changes
# VERY IMPORTTANT Leave a blank line in the start and end of the  cir file

#%% inputs

# dendritic bias current
dIde = 5
Ide_vec = np.arange(210,300+dIde,dIde) # uA

# dendritic integration loop saturation capacity bias current
dIsc = 5
Isc_vec = np.arange(70,100+dIsc,dIsc)

# activity flux input
max_flux = p['Phi0']/2
Ma = np.sqrt(200e-12*10.3e-12)
Ia_max = 1e6*max_flux/Ma

# jtl inductors
Ljtl1 = 20.6 # pH

#%%
rate = WRSpice.WRSpice(tempCirFile = 'dend__4jj__lin_ramp__100uA', OutDatFile = 'dend_4jj_lin_ramp_100uA_') # no file extentions
rate.pathWRSSpice = '/raid/home/local/xictools/wrspice.current/bin/wrspice'  # for running on grumpy
# rate.pathWRSSpice='C:/usr/local/xictools/wrspice/bin/wrspice.bat' # for local computer

rate.save = 'i(L2) v(5)' # list all vectors you want to save
rate.pathCir = ''


# FileFormat = 'Ib{:05.2f}uA_Ip{:05.2f}uA_rmp{:1d}'
FileFormat = 'Ide{:06.2f}uA_Isc{:06.2f}uA_Ljtl1{:5.2f}pH'
# for aa in [-1,1]: # don't want both sides of in activity flux; negative gives undesirable response shape
for ii in range(len(Ide_vec)):
    Ide = Ide_vec[ii]
    for jj in range(len(Isc_vec)):
        Isc = Isc_vec[jj]
        
        print('ii = {} of {}; jj = {} of {}'.format(ii+1,len(Ide_vec),jj+1,len(Isc_vec)))            
            
        FilePrefix = FileFormat.format(Ide,Isc,Ljtl1)
        rate.FilePrefix = FilePrefix
                 
        Params = {          
        'Isc':'{:6.2f}u'.format(np.round(Isc,2)),
        'Ide':'{:6.2f}u'.format(np.round(Ide,2)),
        'Ia':'{:6.2f}u'.format(np.round(Ia_max,2)),
        'Ljtl1':'{:5.2f}p'.format(Ljtl1)
        }
    
        rate.Params = Params
        rate.stepTran = '10p'
        rate.stopTran = '102n'
        rate.doAll()
