import zlib
import pickle
import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse
import scipy.io as sio
import os
import sys
from cxsystem2.core.tools import write_to_file as wtf
from brian2.units import *
import datetime
import pandas as pd
import elephant as el
from neo.core import AnalogSignal
import quantities as pq
import pdb

def parsePath(path,filename, type='results'):
    filename_out = None
    if path is None:
        path = './'
        
    if filename is None:
        filename_out = _newest(path,type=type)
    
    if not filename_out:
        filename_out = os.path.join(path,filename)
    return filename_out

def close():
    plt.close('all')

def figsave(figurename='', myformat='png', suffix=''):
    
    # Saves current figure to working dir

    if myformat[0] == '.':
        myformat=myformat[1:]

    filename, file_extension = os.path.splitext(figurename)

    filename = filename + suffix

    if not file_extension:
         file_extension = '.' + myformat

    if not figurename:
        figurename = 'MyFigure' + file_extension
    else:
        figurename = filename + file_extension

    plt.savefig(figurename, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=file_extension[1:],
        transparent=False, bbox_inches=None, pad_inches=0.1,
        metadata=None)

def getData(filename):

    # If extension is .gz, open pickle, else assume .mat
    filename_root, filename_extension = os.path.splitext(filename)
    if 'gz' in filename_extension:
        fi = open(filename, 'rb')
        try:
            data_pickle = zlib.decompress(fi.read())
            data = pickle.loads(data_pickle)
        except:
            with open(filename, 'rb') as data_pickle:
                data = pickle.load(data_pickle)
    elif 'mat' in filename_extension:
        data = {}
        sio.loadmat(filename,data) 
    else:
        raise TypeError('U r trying to input unknown filetype, aborting...')

    return data

def _getDistance(positions,index_position=None):
    '''Calculates distances between neurons. If index position is given, 
    only distances to this neuron is calculated. Without index position, 
    all distance pairs will be calculated.
    Assumes positions as a list or numpy array of complex coordinates'''
    
    positions_array = np.asarray(positions)
    assert len(positions_array) > 1, 'At least two positions necessary for distance'

    # Calculate distance btw one cell and all other cells
    def dist(index_position,other_positions):
        # Assumes arrays of complex numbers
        d = np.sqrt((np.real(index_position)-np.real(other_positions))**2 + 
            (np.imag(index_position)-np.imag(other_positions))**2)
        return d

    # Check whether index_position exists
    if index_position is not None:
        # Calculate distance between index neuron and other neurons
        distance = dist(index_position,positions_array)
    # If only two positions are given
    elif len(positions_array)==2:
        distance = dist(positions_array[0],positions_array[1])
    # Otherwise
    else:
        # Init result matrix
        distance = np.zeros([len(positions_array),len(positions_array)])

        # Loop index neuron positions
        for idx, index_position in enumerate(positions_array):
            ## Create vector array of other positions
            #other_positions = np.hstack([positions_array[:idx],positions_array[idx+1:]]) 
            # Calculate distance between index neuron and all neurons including itself
            distance[idx,:] = dist(index_position,positions_array)
           
    return distance

def _getNeuronIndex(data, neuron_group, position=0+0j):
    neuron_index=data['positions_all']['w_coord'][neuron_group].index(position)
    return neuron_index

def getLambda(D,MF):
    # Length constant according to Schwabe et al J Neurosci 2006
    # D is diameter of mean axonal length in mm along cx.
    # MF is magnification factor (2.3 mm/deg for V1 at 5 deg ecc; 0.45 mm/deg for MT at 5 deg)
    delta_x = 0.5 * D * MF**-1
    l = -1 * delta_x**-1 * np.log(0.05)
    return l

def _createPositions(distance_between_neurons, cx_radius, ndims=2, coordinate_system='w', zero_first=False):
    
    # Define grid. Put RF center to index 0.
    n_neurons_per_row = int(np.ceil((2 * cx_radius) / distance_between_neurons))
    # Check for even N neurons
    if not np.mod(n_neurons_per_row,2):
         n_neurons_per_row += 1
    center_position = 0+0j # assuming center at 0+0j
    
    if ndims==2:
        # Assuming circular grid, n rows = n columns
        n_neurons_per_column = n_neurons_per_row
        positions_real = np.linspace(-cx_radius,cx_radius,n_neurons_per_row)
        positions_imag = np.linspace(-cx_radius,cx_radius,n_neurons_per_column)
    elif ndims==1:
        positions_real = np.linspace(-cx_radius,cx_radius,n_neurons_per_row)
        #Place all cells to y=0
        positions_imag = 0
    else:
        raise NotImplementedError('Number of dimensions is not 1 or 2')
    
    
    positions_grid = np.meshgrid(positions_real,positions_imag)
    positions_grid[1] = positions_grid[1] * 1j
    positions_grid_np_array = positions_grid[0] + positions_grid[1]
    positions_grid_np_array_flat = positions_grid_np_array.flatten()

    # Cut circle according to cx_radius distance
    distances_to_center = _getDistance(positions_grid_np_array_flat, index_position=center_position)
    positions = positions_grid_np_array.flatten()[distances_to_center<=cx_radius]

    if zero_first:
        # Find center
        center_index_tuple=np.where(positions == center_position)
        assert len(center_index_tuple) == 1, 'Other than 1 values equal center value'
        center_index = int(center_index_tuple[0])
        positions_zero_first = np.hstack((positions[center_index],
                                        positions[:center_index],
                                        positions[center_index+1:]))
        positions=positions_zero_first

    if coordinate_system=='z':
        # magnification factor at 5 deg according to Schwabe 2006
        M=2.3
        positions_z = positions/M
        positions = positions_z
    elif coordinate_system=='w':
        pass
    else:
        raise NotImplementedError('Unknown coordinate system')

    positions_out=[(pos) for pos in positions] 

    # Return positions in complex coordinates
    return positions_out

def buildSchwabePositions(ndims=1,group_keys=None,group_values=None):
    
    # Three first groups are in V1 (161 neurons) and the fourth is in V5 (33)
    if group_keys is None or group_values is None:
        group_keys=['NG0_relay_vpm', 'NG1_SS_L2', 'NG2_MC_L2', 'NG3_SS_L6', 'NG4_L2_SS_autoconn_L2']
        group_values=[(0.23,18.4), (0.23,18.4), (0.23,18.4), (1.15,18.4), (0.23,18.4)]
        print("Using default group names and values")
    
    positions_w = {}
    positions_z = {}
    for group_key, group_value in zip(group_keys, group_values):
        positions_w[group_key] = _createPositions(group_value[0],group_value[1],ndims=ndims, coordinate_system='w') 
        positions_z[group_key] = _createPositions(group_value[0],group_value[1],ndims=ndims, coordinate_system='z') # Dummy, code needs
    coord_dict={'w_coord':positions_w, 'z_coord':positions_z}
    data_positions={'positions_all':coord_dict}
    return data_positions

'''def buildSchwabeConnections(data):

    # Replace id only connections
    connection_strengths = 7e-10
    data['relay_vpm__to__SS_L2_soma']['data'] = connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG1_SS_L2']).size))
    data['SS_L2__to__L2_SS_autoconn_L2_soma']['data'] = connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG1_SS_L2']).size))
    data['L2_SS_autoconn_L2__to__SS_L2_soma']['data'] = connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG1_SS_L2']).size))
    data['MC_L2__to__SS_L2_soma']['data'] = connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG1_SS_L2']).size))
    data['MC_L2__to__MC_L2_soma']['data'] = connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG2_MC_L2']).size))


    Add id to id + lambda connections
    data['SS_L4__to__SS_L4_soma']['data'] = data['SS_L4__to__SS_L4_soma']['data'] + \
                                                connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG1_SS_L4']).size))
    
    This covers both local as well as lateral excitation of inhibitory neurons. They have the same value in Schwabe 2006 model.
    data['SS_L2__to__MC_L2_soma']['data'] = data['SS_L2__to__MC_L2_soma']['data'] + \
                                                connection_strengths * sparse.csr_matrix(
                                                np.identity(np.asarray(data['positions_all']['w_coord']['NG2_MC_L2']).size))

    return data'''

def saveSchwabeData(data, base_filename,dir_name=None):
    filename = base_filename + '.gz'
    if dir_name is None:
        dir_name = os.getcwd()
    filenameWithPath = os.path.join(dir_name, filename)
    wtf(filenameWithPath, data)

def _newest(path='./',type='connections'):
    # type = {'connections','results'}
    files = [f for f in os.listdir(path) if type in f]
    paths = [os.path.join(path, basename) for basename in files]
    filename = max(paths, key=os.path.getctime)
    # fullfile = os.path.join(path,filename)
    # return fullfile
    return filename

def showLatestConnections(path='./',filename=None, hist_from=None, savefigname=''):

    filename = parsePath(path,filename, type='connections')
    
    print(filename)
    data = getData(filename)

    # Visualize
    # Extract connections from data dict
    list_of_connections = [n for n in data.keys() if '__to__' in n]
 
    # Pick histogram data
    if hist_from is None:
        hist_from = list_of_connections[-1]

    print(list_of_connections)
    n_images=len(list_of_connections)
    n_columns = 2
    n_rows = int(np.ceil(n_images/n_columns))
    fig, axs = plt.subplots(n_rows, n_columns)
    axs = axs.flat
    for ax, connection in zip(axs,list_of_connections):
        im = ax.imshow(data[connection]['data'].todense())
        ax.set_title(connection, fontsize=10)
        fig.colorbar(im, ax=ax)
    data4hist = np.squeeze(np.asarray(data[hist_from]['data'].todense().flatten()))
    data4hist_nozeros = np.ma.masked_equal(data4hist,0)
    axs[(n_rows * n_columns)-1].hist(data4hist_nozeros)
    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def showLatestVm(path='./',filename=None, savefigname=''):

    filename = parsePath(path,filename, type='results')

    print(filename)
    data = getData(filename)
    # Visualize
    # Extract connections from data dict
    list_of_results = [n for n in data['vm_all'].keys() if 'NG' in n]

    print(list_of_results)
    n_images=len(list_of_results)
    n_columns = 2
    n_rows = int(np.ceil(n_images/n_columns))

    t=data['vm_all'][list_of_results[0]]['t']
    # time_interval=[2000, 4000]
    time_interval=[0, 2000]

    fig, axs = plt.subplots(n_rows, n_columns)
    axs = axs.flat

    for ax, results in zip(axs,list_of_results):
        N_monitored_neurons = data['vm_all'][results]['vm'].shape[1]
        N_neurons = len(data['positions_all']['w_coord'][results])
        # If all neurons are monitored, show center and it's neighborghs, otherwise, show all.
        if N_monitored_neurons == N_neurons: 
            # neuron_index_center=data['positions_all']['w_coord'][results].index(0+0j)
            neuron_index_center = _getNeuronIndex(data, results, position=0+0j)
            
            im = ax.plot(t[time_interval[0]:time_interval[1]], 
                        data['vm_all'][results]['vm'][time_interval[0]:time_interval[1],
                        neuron_index_center-1:neuron_index_center+2])
        else:
            im = ax.plot(t[time_interval[0]:time_interval[1]], 
                        data['vm_all'][results]['vm'][time_interval[0]:time_interval[1],:])
        ax.set_title(results, fontsize=10)

    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def showLatestG(path='./',filename=None, savefigname=''):

    filename = parsePath(path,filename, type='results')

    print(filename)
    data = getData(filename)
    # Visualize
    # Extract connections from data dict
    list_of_results_ge = [n for n in data['ge_soma_all'].keys() if 'NG' in n]
    list_of_results_gi = [n for n in data['gi_soma_all'].keys() if 'NG' in n]

    print(list_of_results_ge)
    n_images=len(list_of_results_ge)
    n_columns = 2
    n_rows = int(np.ceil(n_images/n_columns))

    t=data['ge_soma_all'][list_of_results_ge[0]]['t']
    time_interval=[2000, 4000]

    fig, axs = plt.subplots(n_rows, n_columns)
    axs = axs.flat

    for ax, results in zip(axs,list_of_results_ge):
        im = ax.plot(t[time_interval[0]:time_interval[1]], 
                     data['ge_soma_all'][results]['ge_soma'][time_interval[0]:time_interval[1],0:-1:20])
        ax.set_title(results + ' ge', fontsize=10)

    fig2, axs2 = plt.subplots(n_rows, n_columns)
    axs2 = axs2.flat

    for ax2, results2 in zip(axs2,list_of_results_gi):
        im = ax2.plot(t[time_interval[0]:time_interval[1]], 
                     data['gi_soma_all'][results2]['gi_soma'][time_interval[0]:time_interval[1],0:-1:20])
        ax2.set_title(results2 + ' gi', fontsize=10)

    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def showLatestI(path='./',filename=None, savefigname=''):

    filename = parsePath(path,filename, type='results')

    print(filename)
    data = getData(filename)
    # Visualize
    # Extract connections from data dict
    list_of_results_ge = [n for n in data['ge_soma_all'].keys() if 'NG' in n]
    list_of_results_gi = [n for n in data['gi_soma_all'].keys() if 'NG' in n]
    list_of_results_vm = [n for n in data['vm_all'].keys() if 'NG' in n]

    El = -65 * mV
    gl = 50 * nS
    Ee = 0 * mV
    Ei = -75 * mV

    print(list_of_results_ge)
    print(f'Assuming El = {El:6.4f} * mV, gl = {gl:6.4f} * nS, Ee = {Ee:6.4f} * mV, Ei = {Ei:6.4f} * mV\nDig from physiology df if u start messing with neuron types')

    n_images=len(list_of_results_ge)
    n_columns = 2
    n_rows = int(np.ceil(n_images/n_columns))

    t=data['ge_soma_all'][list_of_results_ge[0]]['t']
    time_interval=[2000, 4000]

    fig, axs = plt.subplots(n_rows, n_columns)
    axs = axs.flat

    for ax, results_ge, results_gi, results_vm in zip(  axs,list_of_results_ge, \
                                                        list_of_results_gi,list_of_results_vm):

        N_monitored_neurons = data['vm_all'][results_vm]['vm'].shape[1]
        N_neurons = len(data['positions_all']['w_coord'][results_vm])

        ge= data['ge_soma_all'][results_ge]['ge_soma']
        gi= data['gi_soma_all'][results_gi]['gi_soma']
        vm= data['vm_all'][results_vm]['vm']
        I_total = gl * (El - vm) + ge * (Ee - vm) + gi * (Ei - vm)

        if N_monitored_neurons == N_neurons: 
            # neuron_index_center=data['positions_all']['w_coord'][results_vm].index(0+0j)
            neuron_index_center = _getNeuronIndex(data, results_vm, position=0+0j)
            ax.plot(t[time_interval[0]:time_interval[1]], 
                        I_total[time_interval[0]:time_interval[1],
                        neuron_index_center])
        else:
            ax.plot(t[time_interval[0]:time_interval[1]], 
                        I_total[time_interval[0]:time_interval[1],:])

        ax.set_title(results_vm + ' I', fontsize=10)

        I_total_mean = np.mean(I_total[time_interval[0]:time_interval[1]] / namp)
        I_total_mean_str = f'mean I = {I_total_mean:6.2f} nAmp'
        ax.text(0.05, 0.95, I_total_mean_str, fontsize=10, verticalalignment='top', transform=ax.transAxes)

    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def _getI(neuron_group,data):

    # Extract connections from data dict
    # list_of_results_ge = [n for n in data['ge_soma_all'].keys() if 'NG' in n]
    # list_of_results_gi = [n for n in data['gi_soma_all'].keys() if 'NG' in n]
    # list_of_results_vm = [n for n in data['vm_all'].keys() if 'NG' in n]
    list_of_results_ge = [n for n in data['ge_soma_all'].keys() if neuron_group in n]
    list_of_results_gi = [n for n in data['gi_soma_all'].keys() if neuron_group in n]
    list_of_results_vm = [n for n in data['vm_all'].keys() if neuron_group in n]
    assert list_of_results_ge == list_of_results_gi == list_of_results_vm, "different N neuron groups monitored, can not calculate current, aborting..."

    if len(list_of_results_ge)==0:
        #No currents monitored for thir group
        I_total_mean = 0 # to avoid stupid error later
        I_excitatory_mean = 0
        I_inhibitory_mean = 0
        return I_total_mean, I_excitatory_mean, I_inhibitory_mean
    # These data are available in physiology_configuration file (inside results datafile), but hard coded here for the time being
    El = -65 * mV
    gl = 50 * nS
    Ee = 0 * mV
    Ei = -75 * mV

    time_interval=[2000, 4000]

    for results_ge, results_gi, results_vm in zip(  list_of_results_ge, \
                                                        list_of_results_gi,list_of_results_vm):

        N_monitored_neurons = data['vm_all'][results_vm]['vm'].shape[1]
        N_neurons = len(data['positions_all']['w_coord'][results_vm])

        ge= data['ge_soma_all'][results_ge]['ge_soma']
        gi= data['gi_soma_all'][results_gi]['gi_soma']
        vm= data['vm_all'][results_vm]['vm']
        I_total = gl * (El - vm) + ge * (Ee - vm) + gi * (Ei - vm)
        I_excitatory = ge * (Ee - vm)
        I_inhibitory = gi * (Ei - vm)

        if N_monitored_neurons == N_neurons: 
            neuron_index_center = _getNeuronIndex(data, results_vm, position=0+0j)
        I_total_mean = np.mean(I_total[time_interval[0]:time_interval[1]] / namp)
        I_excitatory_mean = np.mean(I_excitatory[time_interval[0]:time_interval[1]] / namp)
        I_inhibitory_mean = np.mean(I_inhibitory[time_interval[0]:time_interval[1]] / namp)
        
    return I_total_mean, I_excitatory_mean, I_inhibitory_mean

def showLatestSpatial(path='./',filename=None,sum_length=1, savefigname=''):

    filename = parsePath(path,filename, type='results')

    print(filename)
    data = getData(filename)

    # Visualize
    coords='w_coord'
    # Extract connections from data dict
    list_of_results = [n for n in data['spikes_all'].keys() if 'NG' in n]

    print(list_of_results)
    n_images=len(list_of_results)
    n_columns = 2
    ylims=np.array([-18.4,18.4])
    n_rows = int(np.ceil(n_images/n_columns))

    width_ratios = np.array([3,1,3,1])
    fig, axs = plt.subplots(n_rows, n_columns * 2, gridspec_kw={'width_ratios': width_ratios})
    axs = axs.flat
    # flat_list_of_results = [item for sublist in list_of_results for item in sublist]

    # Spikes
    for ax1, results in zip(axs[0:-1:2],list_of_results):
        # im = ax1.scatter(data['spikes_all'][results]['t'], data['spikes_all'][results]['i'],s=1)
        position_idxs=data['spikes_all'][results]['i']
        im = ax1.scatter(   data['spikes_all'][results]['t'], 
                            np.real(data['positions_all'][coords][results])[position_idxs],
                            s=1)
        ax1.set_ylim(ylims)
        ax1.set_title(results, fontsize=10)
    # Summary histograms
    for ax2, results in zip(axs[1:-1:2],list_of_results):
        # position_idxs=data['spikes_all'][results]['i']
        #Rearrange and sum neurons in sum_length bin width
        if sum_length is None:
            sum_length = 3 # How many neurons to sum together
        full_vector=data['spikes_all'][results]['count'].astype('float64')
        necessary_length = int(np.ceil(len(full_vector)/sum_length)) * sum_length
        n_zeros = necessary_length - len(full_vector)
        full_vector_padded = np.pad(full_vector, (0, n_zeros), 'constant', constant_values=(np.NaN, np.NaN))
        rearranged_data = np.reshape(full_vector_padded,[sum_length, int(necessary_length/sum_length) ], order='F')
        summed_rearranged_data = np.nansum(rearranged_data, axis=0)
        full_vector_pos = np.real(data['positions_all'][coords][results])
        full_vector_pos_padded = np.pad(full_vector_pos, (0, n_zeros), 'constant', constant_values=(np.NaN, np.NaN))
        rearranged_pos = np.reshape(full_vector_pos_padded,[sum_length, int(necessary_length/sum_length) ], order='F')
        mean_rearranged_pos = np. nanmean(rearranged_pos, axis=0)

        firing_frequency = summed_rearranged_data / (sum_length * data['runtime'])

        pos = mean_rearranged_pos
        im = ax2.barh(pos, firing_frequency, height=1.0)
        ax2.set_ylim(ylims)
        # ax2.axes.get_yaxis().set_visible(False)

    # Shut down last pair if uneven n images
    if np.mod(n_images,2):
        axs[-2].axis('off')
        axs[-1].axis('off')
    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def showLatestSpectra(path='./',filename=None, savefigname=''):

    filename = parsePath(path,filename, type='results')

    print(filename)
    data = getData(filename)

    # Visualize
    coords='w_coord'
    # Extract connections from data dict
    list_of_results = [n for n in data['spikes_all'].keys() if 'NG' in n]

    print(list_of_results)
    n_images=len(list_of_results)
    n_columns = 2
    nbins=400
    duration = data['runtime']
    bar_width = duration / nbins
    freqCutoff = 200 # Hz
    ylims=np.array([-18.4,18.4])
    n_rows = int(np.ceil(n_images/n_columns))

    width_ratios = np.array([2,2])

    fig, axs = plt.subplots(n_rows * 2, n_columns, gridspec_kw={'width_ratios': width_ratios})
    axs = axs.flat

    for ax1, results in zip(axs[0:-1:2],list_of_results):
        # im = ax1.scatter(data['spikes_all'][results]['t'], data['spikes_all'][results]['i'],s=1)
        counts, bins = np.histogram(data['spikes_all'][results]['t']/second, range=[0,duration], bins=nbins)
        Nneurons = data['number_of_neurons'][results]
        firing_rates = counts * (nbins / (duration * Nneurons))
        im = ax1.bar(bins[:-1], firing_rates, width=bar_width)
        ax1.set_title(results, fontsize=10)

    for ax2, results in zip(axs[1:-1:2],list_of_results):
        try:
            Vms=data['vm_all'][results]['vm']
        except:
            continue
        sigarr = AnalogSignal(Vms, units='mV', sampling_rate=10000*pq.Hz)
        freqs, psd = el.spectral.welch_psd(sigarr, n_segments=8, len_segment=None, frequency_resolution=None, fs=10000, nfft=None, scaling='density', axis=- 1)
        
        # Cut to desired freq range
        freqs_for_plotting = freqs[freqs < freqCutoff]
        psd_for_plotting = psd[:,freqs < freqCutoff]

        # Calculate mean psd
        psd_for_plotting_mean = np.mean(psd_for_plotting, axis=0)
 
        ax2.plot(freqs_for_plotting, psd_for_plotting.T, color='gray', linewidth=.25)
        ax2.plot(freqs_for_plotting, psd_for_plotting_mean, color='black', linewidth=1)

    # Shut down last pair if uneven n images
    if np.mod(n_images,2):
        axs[-2].axis('off')
        axs[-1].axis('off')
    if savefigname:
        figsave(figurename=savefigname)

    plt.show()

def showLatestASF(path='./',timestamp=None,sum_length=1, fixed_y_scale=True, data_type='spikes', savefigname=''):
    '''
    Show ASF curves for single value or for an array search across one independent variable.
    At the moment, no more than 1000 values are adviced for one search,
    e.g. with 20 trials per run and 10 ASF sized this means 5 values for the array of independent variable.
    Dimensions are (ASF_size, search_variable, trial)
    '''
    
    assert data_type=='spikes' or data_type=='current', "Unknown data type, should be 'spikes' or 'current', aborting"
    
    if timestamp is None:
        filename = parsePath(path,filename=None, type='results')
        today = str(datetime.date.today()).replace('-','')
        start_index = filename.find(today)
        end_index = start_index + 16
        timestamp = filename[start_index:end_index]

    # Get all files with the same timestamp
    # list all files with the timestamp
    all_files = os.listdir(path)
    files_correct_timestamp = [files for files in all_files if timestamp in files]
    # skip metadata and connections files, get other filenames
    result_files_for_ASF_step1 = [files for files in files_correct_timestamp if 'metadata' not in files]
    result_files_for_ASF = [files for files in result_files_for_ASF_step1 if 'connections' not in files]

    # Sort the result files to increasing center size
    metadata_file = [files for files in files_correct_timestamp if 'metadata' in files]
    assert len(metadata_file)<=1, "Multiple metadatafiles, don't know what to do, aborting"
    assert len(metadata_file)==1, "No metadatafile, cannot do ASF from single file, or from different runs"
    metadata_df=getData(os.path.join(path, metadata_file[0]))

    # Test if multiple trials per run and one or two dimensions (first is always ASF size)
    # trials_per_config = Number of files / Number of unique parameters
    assert 'Dimension-3 Value' not in metadata_df.keys(), 'Cannot handle more than 2 dims for one ASF analysis, aborting'
    ASF_array_length = metadata_df['Dimension-1 Value'].unique().size
    if 'Dimension-2 Value' in metadata_df.keys():
        print('\nTwo-dimensional array run')
        search_variable_name = metadata_df['Dimension-2 Parameter'].unique()[0]
        search_variable_array = metadata_df['Dimension-2 Value'].unique()
        search_variable_array_length = search_variable_array.size
        trials_per_config =  int(   metadata_df['Full path'].size /    
                                    (ASF_array_length * search_variable_array_length))
    else:
        search_variable_name = 'no'
        search_variable_array = ['--']
        search_variable_array_length = 1
        trials_per_config =  int(metadata_df['Full path'].size / metadata_df['Dimension-1 Value'].unique().size)

    # assert 'act' in result_files_for_ASF[0], 'I was expecting "act" in filename. This might not be ASF data, aborting'
    
    contrast_types = np.array(['ON', 'OFF', 'ACT'])
    contrast_type_idx = [contrast_type in result_files_for_ASF[0] for contrast_type in contrast_types]
    contrast_type = contrast_types[np.array(contrast_type_idx)][0]
    assert contrast_type, 'I was expecting "on", "off" or "act" in filename. This might not be ASF data, aborting'


    #The sizes are extracted from filenames, so BE CAREFUL.
    # Assuming annulus data if filename begins with 'ANN'
    if result_files_for_ASF[0][:3]=='ANN':
        assert  result_files_for_ASF[0].count(contrast_type)==3, '''Annulus data should have "act", "on" or "off" three times in results data filename ie 
                                                                    in the spike_times variable, aborting'''
        filename_dict = {}
        for this_file in result_files_for_ASF:
            relevant_substring = this_file[this_file.index(contrast_type):this_file.index('.gz')]
            try:
                active_ranges = eval("np.array([[" + relevant_substring[3:].replace(contrast_type,'],[').replace('-',',') + "]])")
            except:
                relevant_substring = relevant_substring[:relevant_substring.index('_')] # 2-dim run puts _nextDimParmas to end
                active_ranges = eval("np.array([[" + relevant_substring[3:].replace(contrast_type,'],[').replace('-',',') + "]])")

            annulus_distance = (np.mean(active_ranges[0,:]) - active_ranges[1,1] + active_ranges[2,0] - np.mean(active_ranges[0,:])) / 2
            filename_dict[this_file] = annulus_distance

        filename_array_sorted = metadata_df['Full path'].iloc[::-1].values
    
    # This should be ASF data
    elif result_files_for_ASF[0][:3]=='ASF':
        if trials_per_config==1:
            # Get ASF radius (SIC) from filenames: pick str btw 'act' and '-1', eval the difference, multiply by -1 to get it positive
            filename_dict = {k:(-1 * eval(k[k.find(contrast_type) + len(contrast_type) : k.find('-1.gz')]))/2 for k in result_files_for_ASF}
        elif trials_per_config>1:
            filename_dict = {k:(-1 * eval(k[k.find(contrast_type) + len(contrast_type) : k.find('-1_')]))/2 for k in result_files_for_ASF}
        filename_array_sorted = metadata_df['Full path'].values

    else:
        assert 0, "Neither ASF nor ANN in filename, aborting"


    # Get unique ASF_size values with set command
    ASF_x_axis_values = sorted(set(filename_dict.values()))

    coords='w_coord'

    # Get neuron group names and init ASF_dicts
    data = getData(os.path.join(path, result_files_for_ASF[0]))
    list_of_results = [n for n in data['spikes_all'].keys() if 'NG' in n]
    if data_type == 'spikes':
        ASF_dict = {k:np.zeros([ASF_array_length,search_variable_array_length, trials_per_config]) for k in list_of_results}
    elif data_type == 'current':
        ASF_dict = {k:np.zeros([ASF_array_length,search_variable_array_length, trials_per_config, 3]) for k in list_of_results}

    # Loop for data
    trial = 0
    search_variable = 0
    ASF_size = 0

    epoch_duration = np.max(data['time_vector']) - np.min(data['time_vector'])

    # for filename in filename_dict_sorted:
    for filename in filename_array_sorted:
        data = getData(filename)
        # print(f'I am happy for data {filename}')
        # For each neuron group, get spike frequencies for neurons of interest, accept sum_length
        for neuron_group in list_of_results:
            # Pick center idx
            center_index = _getNeuronIndex(data, neuron_group, position=0+0j)
            neuron_indices = np.arange(center_index - np.floor(sum_length/2), center_index + np.ceil(sum_length/2)).astype('int')
            
            # Select data, add to nd array
            if data_type == 'spikes':
                full_vector=data['spikes_all'][neuron_group]['count'].astype('float64')
                firing_frequency = np.mean(full_vector[neuron_indices]) / epoch_duration
                ASF_dict[neuron_group][ASF_size, search_variable, trial] = firing_frequency
            elif data_type == 'current':
                I_total_mean, I_excitatory_mean, I_inhibitory_mean = _getI(neuron_group, data)
                ASF_dict[neuron_group][ASF_size, search_variable, trial, 0] = I_total_mean
                ASF_dict[neuron_group][ASF_size, search_variable, trial, 1] = I_excitatory_mean
                ASF_dict[neuron_group][ASF_size, search_variable, trial, 2] = I_inhibitory_mean *-1 ## Inverting inhibitory currents to positive


        trial += 1
        if trial == trials_per_config: 
            trial = 0; search_variable += 1
        if search_variable == search_variable_array_length: 
            search_variable = 0; ASF_size += 1

    # For ANN data the ANN size is inverted above, resulting in inverting search variable and trial, too.
    # Here, the search variable will be inverted. I also invert the trial number for consistency.
    if result_files_for_ASF[0][:3]=='ANN':
        ASF_dict_inv = {}
        for neuron_group in list_of_results:
            ASF_dict_inv[neuron_group] = np.flip(ASF_dict[neuron_group], axis=(1,2)) # dim 1 = search variable, dim 2 = trial
        ASF_dict = ASF_dict_inv

    # Visualize
    n_images=len(list_of_results)
    n_columns = 2
    n_rows = int(np.ceil(n_images/n_columns))
    
    # Enable interactive mode, ie don't block ipython with plt.show()
    plt.ion() 
    if data_type == 'spikes':
        labels = ('firing_frequency')
    elif data_type == 'current':
        labels = ('I_total_mean', 'I_excitatory_mean','I_inhibitory_mean')

    # ASF is defined for center receptive field. I leave here the option to show the
    # spatial dimension in the second subplot, including eg the parameters of DoG fit
    # along space. This gives us reflection of the FB in case the HO model is not perfect.
    # width_ratios = np.array([3,1,3,1])
    width_ratios = np.array([3,3])
    if fixed_y_scale:
        y_scales=dict(zip(list_of_results, [None]*len(list_of_results)))
        for neuron_group in list_of_results:
            y_max = np.max(ASF_dict[neuron_group])
            # Get y_scales with 10% over max not to obliterate high values in the plot
            y_scales[neuron_group] = y_max + np.ceil(y_max /  10) 

    # Loop for 2nd dim, ie for the primary independent search variable; the ASF is fixed
    # One figure for one search variable value
    for search_variable_index in np.arange(search_variable_array_length):
        # fig, axs = plt.subplots(n_rows, n_columns * 2, gridspec_kw={'width_ratios': width_ratios}, num=search_variable_index)
        fig, axs = plt.subplots(n_rows, n_columns, gridspec_kw={'width_ratios': width_ratios}, num=search_variable_index)
        fig.suptitle(f'Searching for {search_variable_name} parameter, value {search_variable_array[search_variable_index]}', fontsize=16)
        axs = axs.flat

        # Spikes
        # for ax1, neuron_group in zip(axs[0:-1:2],list_of_results):
        for ax1, neuron_group in zip(axs,list_of_results):
            # pick data
            data_to_plot = np.squeeze(ASF_dict[neuron_group][:,search_variable_index,:])
            if data_to_plot.ndim == 1:
                data_to_plot_mean = data_to_plot
            else:
                data_to_plot_mean = np.squeeze(data_to_plot.mean(axis=1))

            # plot
            if data_to_plot_mean.ndim == 1:
                ax1.plot(ASF_x_axis_values,data_to_plot, 'k', linewidth=0.5, alpha=0.1)
                ax1.plot(ASF_x_axis_values,data_to_plot_mean, 'k', linewidth=2, label=labels)
            else:
                ax1.plot(ASF_x_axis_values,data_to_plot_mean, linewidth=2, label=labels)

            if fixed_y_scale:
                ax1.set_ylim(bottom=0, top=y_scales[neuron_group])
            ax1.set_title(neuron_group, fontsize=10)

        # Shut down last pair if uneven n images
        if np.mod(n_images,2):
            # axs[-2].axis('off')
            axs[-1].axis('off')

        handles, foo = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower right')
        if savefigname:
            figsave(figurename=savefigname, suffix=f'{search_variable_array[search_variable_index]}')

    plt.show()

def createASFset(start_stim_radius=0.1,end_stim_radius=8,units='deg',Nsteps=5,show_positions=True):
    
    M=2.3 # M factor of macaque V1 at 5 deg ecc

    # Turn units to mm cortex
    if units == 'deg':
        start_stim_radius = start_stim_radius * M
        end_stim_radius = end_stim_radius * M
    else:
        assert units == 'mm', "Unknown units, aborting"
    
    # Name neuron group_keys
    group_keys=['NG0_relay_vpm', 'NG1_SS_L2', 'NG2_MC_L2', 'NG3_SS_L6', 'NG4_L2_SS_autoconn_L2']

    # Set outer radius of cells in mm (V1 widest, V5 somewhat smaller to see the edge effect)
    V1_outer_radius = 18.4
    V1_distance_btw_neurons = 0.23
    V5_outer_radius = 18.4
    V5_distance_btw_neurons = 1.15

    # Create positions (start, end, step) with buildSchwabePosisions
    input_radii = np.round( np.linspace(start_stim_radius, end_stim_radius, Nsteps),
                            decimals=1)

    # Create dictionary holding all input radii
    ASF_positions_dict={}
    for input_radius in input_radii:
        # Calculate rounded input radii for same positions as in V1
        input_radius_match_V1 = np.ceil(input_radius / V1_distance_btw_neurons) * V1_distance_btw_neurons
        group_values=[  (V1_distance_btw_neurons,input_radius_match_V1), 
                        (V1_distance_btw_neurons,V1_outer_radius), 
                        (V1_distance_btw_neurons,V1_outer_radius), 
                        (V5_distance_btw_neurons,V5_outer_radius), 
                        (V1_distance_btw_neurons,V1_outer_radius)]

        data_positions = buildSchwabePositions( ndims=1,group_keys=group_keys,
                                                group_values=group_values)

        key='in_radius_' + str(input_radius)
        ASF_positions_dict[key]=data_positions

        if show_positions:
            plt.figure()
            coords='z_coord' # ['w_coord' | 'z_coord']
            counter = 0
            for neuron_group in group_keys:
                points = data_positions['positions_all'][coords][neuron_group]
                plt.scatter(np.real(points),counter + np.imag(points),s=1,)
                plt.grid(color='k', linestyle='--',axis='x')
                counter += 1
            plt.show()

        # Save position files with compact names by calling saveSchwabeData
        saveSchwabeData(data_positions, key, dir_name='../connections')

def pp_df_memory(df):
    BYTES_TO_MB_DIV = 0.000001
    mem = round(df.memory_usage().sum() * BYTES_TO_MB_DIV, 3) 
    print("Memory usage is " + str(mem) + " MB")

def pp_df_full(df):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None): 
        print(df)

def pp_obj(obj):
    from IPython.lib.pretty import pprint 
    pprint(obj)
    print(f'\nObject size is {sys.getsizeof(obj)} bytes')

def read_neo_block(filename):
    import neo
    reader = neo.io.NixIO(filename=filename, mode='ro')
    block = reader.read_block() 
    return block

if __name__=='__main__':
    path = r'C:\Users\Simo\Laskenta\Models\Grossberg\CxPytestWorkspace\all_results'
    os.chdir(path)
    showLatestSpectra(filename='results_low_500NG4-3_results_20200826_1153084.gz', savefigname='results_low_500NG4-3.eps')