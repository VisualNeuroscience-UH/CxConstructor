'''
Build images of receptive fields and internal representations of simulated data.
'''

import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import os
import sys
from utilities import getData, figsave
import pdb
    
class InternalModel:
    if sys.platform is 'linux':
        root_path = r'/opt2/Laskenta_ssd/Models/SchwabeModel/ASF_cont_bio'
    elif sys.platform.startswith('win'):
        root_path = r'C:\Users\Simo\Laskenta\Models\SchwabeModel\data_pikkuveljeltä\ASF_ConSearch'

    filename = 'ASF_cont_bio_results_20200617_1352495_spike_times100HON65-96-1_custom_weight24nS_000.gz'
    ng_name = 'NG1_SS_L2'
    input_group_prefix = 'NG0'
    myformat = 'eps'

    def __init__(self, ng_name='',filename='',show_neuron_idx=None, savefigname=None):

        if not ng_name:
            ng_name = self.ng_name
        if not filename:
            filename = self.filename

        self.show_neuron_idx = show_neuron_idx
        self.savefigname = savefigname
        self.get_internal_image(ng_name,filename,self.root_path)
    
    def _get_group_names_from_suffix(self, suffixes, all_group_names):
        group_names_list = []
        for suffix in suffixes:
            for group in all_group_names:
                if group[group.find('_') + 1:] in suffix:
                    group_names_list.append(group)

        return group_names_list

    def analyze_groups(self, all_group_names, connection_data):

        input_group_name = [group for group in all_group_names if group.startswith(self.input_group_prefix)][0]
        
        # Find first order groups
        input_name_prefix = input_group_name[input_group_name.find('_') +1 :] 
        # First_order_connections
        connection_names = connection_data.keys()
        first_order_connections = [conn for conn in connection_names if conn.startswith(input_name_prefix)]

        first_order_group_suffix = [conn[conn.find('__to__') + 6:] for conn in first_order_connections]

        first_order_groups_list = self._get_group_names_from_suffix(first_order_group_suffix, all_group_names)
        
        first_order_name_prefix = [ group[group.find('_') +1 :] for group in first_order_groups_list]

        # build dict of second order connections
        all_second_order_connections_dict = {}
        for this_prefix in first_order_name_prefix:
            second_order_connections4this_prefix = [conn for conn in connection_names if conn.startswith(this_prefix)]
            all_second_order_connections_dict[this_prefix]=second_order_connections4this_prefix

        second_order_suffix = []
        for this_first_order_group in all_second_order_connections_dict.keys():
            for this_connection in all_second_order_connections_dict[this_first_order_group]:
                if this_connection.count(this_first_order_group) > 1:
                    continue
                elif 'autoconn' in this_connection:
                    continue
                else:
                    suffix_name =this_connection[this_connection.find('__to__') + 6:]
                    second_order_suffix.append(suffix_name)

        second_order_groups_list = self._get_group_names_from_suffix(second_order_suffix, all_group_names)

        second_order_connections_dict = {}
        # Get incoming connections for each second order group
        for s_group in second_order_groups_list:
            for f_group in first_order_groups_list:

                candidate_connection = f_group[f_group.find('_') + 1 :] + '__to__' + s_group[s_group.find('_') + 1 :]

                valid_connections=[conn for conn in connection_names if conn.startswith(candidate_connection)]
                second_order_connections_dict[s_group] = valid_connections

        group_analysis_dict = {
            'input' : input_group_name,
            'first_order' : first_order_groups_list,
            'second_order' : second_order_groups_list,
            'first_order_connections' : first_order_connections,
            'second_order_connections_dict' : second_order_connections_dict
            }
        return group_analysis_dict

    def show_rf(self, receptive_fields, ng_name, show_neuron_idx=None, internal_image=None):
            
        input_space = receptive_fields['input_space']

        y, x = np.imag(input_space), np.real(input_space)
        if np.ptp(x) == 0:
            xmin = 0; xmax = 1
        else:
            xmin = x.min(); xmax = x.max()
        if np.ptp(y) == 0:
            ymin = 0; ymax = 1
        else:
            ymin = y.min(); ymax = y.max()     

        if show_neuron_idx:
            conn_scaled = receptive_fields[ng_name]
            z = conn_scaled[show_neuron_idx,:].todense()
            plt.title(f'Receptive field for neuron {show_neuron_idx} in group {ng_name}') 

        else:
            z = internal_image
            plt.title(f'Internal image for group {ng_name}') 


        z_max = np.max(z)

        h = plt.imshow(z.squeeze(), cmap ='Greys', vmin = 0, vmax = z_max, 
                        extent =[xmin, xmax, ymin, ymax], 
                            interpolation ='nearest', origin='lower') 
        plt.colorbar(h) 

        if self.savefigname:
            figsave(figurename=self.savefigname, myformat=self.myformat, suffix=f'_internal_{ng_name}')

        plt.show() 

    def get_rf(self, ng_name, all_group_names, simulation_filename):

        # From connections, get neuron positions and weights from input group to target group
        # In case of second order neurons, get the first-level RF first, then use them to get the
        # second level weights. w_coord is cortical coordinates, z_coord is visual field coordinates.

        connection_filename = simulation_filename.replace('results','connections')
        connection_data = getData(connection_filename)

        group_analysis_dict = self.analyze_groups(all_group_names, connection_data)

        positions_input = connection_data['positions_all']['z_coord'][group_analysis_dict['input']]
        input_space = np.asarray(positions_input)

        receptive_fields = {}
        receptive_fields['input_space'] = input_space

        # Get first order connection rf:s
        for group, connection in zip(group_analysis_dict['first_order'], group_analysis_dict['first_order_connections']):

            conn_mx = connection_data[connection]['data']
            conn_mx_scaled = conn_mx / np.max(conn_mx)
            receptive_fields[group] = np.transpose(conn_mx_scaled)

        # Get second order connection rf:s
        for group in group_analysis_dict['second_order']:
            connections = group_analysis_dict['second_order_connections_dict'][group]
            
            rf2 = np.zeros([1, 1])
            for connection in connections:
                # Matrices are M1 = G1 x IN and M2 = G1 x G2. The correct calculation
                # is M2.T * M1, after which you get M3 = G2 x IN dimension

                G1_name_stub = connection[:connection.find('__to__')]
                G1_name = [g for g in receptive_fields.keys() if g.endswith(G1_name_stub)]
                M1 = receptive_fields[G1_name[0]]
                
                # Scale second order connection
                conn_mx = connection_data[connection]['data']
                M2 = conn_mx / np.max(conn_mx)

                M3 = np.dot(np.transpose(M2), M1)

                try: 
                    rf2 += M3
                except UnboundLocalError: 
                    rf2 = M3 

            # multiply with first order connections RF 2nd = sigma (Conn 2nd * RF 1st)
            receptive_fields[group] = csr_matrix(rf2 / np.max(rf2))

        if self.show_neuron_idx:
            show_neuron_idx = self.show_neuron_idx
            self.show_rf(receptive_fields, ng_name, show_neuron_idx)
        
        return receptive_fields


    def get_internal_image(self, ng_name,simulation_filename,path):
        # Pick neuron group, get rf for each neuron, weigh the rf with response rate from simulation. 
        # Show internal image
        full_filepath = os.path.join(path, simulation_filename)
        data = getData(full_filepath)
        simulation_time = data['runtime']
        all_group_names = data['number_of_neurons'].keys()

        assert ng_name in  all_group_names, 'Neuron group not found'

        # RF matrices are G x IN
        receptive_fields = self.get_rf(ng_name, all_group_names, full_filepath)

        # laske sisäiset edustukset alla kertomalla rf:t aktiivisuuksilla
        ng_spike_freq_data = data['spikes_all'][ng_name]['count'] / simulation_time
        ng_receptive_fields = receptive_fields[ng_name]
        ng_spikes = np.dot(ng_spike_freq_data, ng_receptive_fields.todense())
        internal_image = ng_spikes / np.max(ng_spikes)

        self.show_rf(receptive_fields, ng_name, internal_image=internal_image)
        

if __name__=='__main__':
    mod = InternalModel(ng_name='NG3_SS_L6', show_neuron_idx=None ) 