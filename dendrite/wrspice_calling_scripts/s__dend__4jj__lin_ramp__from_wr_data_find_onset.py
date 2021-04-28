#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

# from soen_sim import input_signal, synapse, dendrite, neuron
from _functions import save_session_data, read_wr_data
from _util import physical_constants

p = physical_constants()
Phi0 = p['Phi0__pH_ns']

plt.close('all')

#%% inputs

# dendritic bias current
dIde = 5
Ide_vec = np.arange(225,300+dIde,dIde) # uA

# dendritic integration loop saturation capacity bias current
dIsc = 5
Isc_vec = np.arange(70,100+dIsc,dIsc)

# activity flux input
max_flux = p['Phi0']/2
Ma = np.sqrt(200e-12*10.3e-12)
Ia_max = 1e6*max_flux/Ma

# jtl inductors
Ljtl1 = 20.6 # pH
num_I_de = len(Ide_vec)

min_peak_height = 182e-6 # units of volts for WR
min_peak_distance = 10 # units of samples

Phi_a_on = np.zeros([len(Ide_vec),len(Isc_vec)])
I_drive_on = np.zeros([len(Ide_vec),len(Isc_vec)])
for ii in range(len(Ide_vec)):
    Ide = Ide_vec[ii]
    for jj in range(len(Isc_vec)):
        Isc = Isc_vec[jj]
        
        print('ii = {} of {}; jj = {} of {}'.format(ii+1,len(Ide_vec),jj+1,len(Isc_vec)))                

        directory = 'wrspice_data/4jj'
        file_name = 'dend_4jj_lin_ramp_100uA_Ide{:06.2f}uA_Isc{:06.2f}uA_Ljtl1{:5.2f}pH.dat'.format(Ide,Isc,Ljtl1)
        j_di_str = 'v(5)'
        # j_di_phase_str = 'v(12)'
        I_a_str = 'L2#branch'
        data_dict = read_wr_data(directory+'/'+file_name)
                    
        # assign data
        time_vec = data_dict['time']
        j_di = data_dict[j_di_str]
        I_a = 1e6*data_dict[I_a_str]
                    
        # find peaks    
        j_di_peaks, _ = find_peaks(j_di, height = min_peak_height, distance = min_peak_distance)
                    
        I_drive_on[ii,jj] = I_a[j_di_peaks[0]]
        Phi_a_on[ii,jj] = Ma*I_a[j_di_peaks[0]]
    
        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(time_vec*1e9,j_di*1e6)
        ax.plot(time_vec[j_di_peaks[0]]*1e9,j_di[j_di_peaks[0]]*1e6,'x')
        ax.set_xlim([0,time_vec[j_di_peaks[1]]*1e9])
        ax.set_xlabel(r'time [ns]')
        ax.set_ylabel(r'$J_{di}$ [$\mu$V]')
        plt.show()
            
#%% save data
save_string = 'dend_4jj_100uA_flux_onset'
data_array = dict()
data_array['Ide_vec'] = Ide_vec # uA
data_array['Isc_vec'] = Isc_vec # uA
data_array['Phi_a_on'] = Phi_a_on # pH ns
save_session_data(data_array,save_string+'.soen',False)

np.set_printoptions(precision=9)
_ts = 'I_drive_on__vec = ['
for ii in range(len(I_drive_on)):
    _ts = '{}{:7.4f},'.format(_ts,I_drive_on[ii])
print('{}]'.format(_ts[0:-1]))
print('Phi_a_on = {}'.format(Phi_a_on))
print('Phi_a_on/Phi0 = {}'.format(Phi_a_on/Phi0))
