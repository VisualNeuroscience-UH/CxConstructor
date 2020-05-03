import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd

from cxsystem2.core.tools import  read_config_file

import pdb

root_path = r'C:\Users\Simo\Laskenta\PythonUtilities\MacV1Buildup'
path_to_tables = os.path.join(root_path, 'tables')
path_to_ni_csv = os.path.join(root_path, 'ni_csv_copy')
path_to_config_files = os.path.join(root_path, 'config_files')


class Config:
    '''
    This class instance contains general objects, concerning all areas. It is here for taking data around and including general methods.
    '''
    def __init__(self, replace_existing_cell_groups=True, anatomy_config_df=None, physiology_config_df=None):
        
        self.replace_existing_cell_groups = replace_existing_cell_groups
        self.anatomy_config_df = anatomy_config_df
        self.physiology_config_df = physiology_config_df

        # Read data from Table 1 (surface area in mm2 for V1, V2 and V5) and Table 2 (N neurons and synapses for V1 in distinct layers) 
        table1_df = self.read_data_from_tables(path_to_tables, 'table1_data.csv')
        self.table1_df = table1_df.set_index('stat') 
        table2_df = self.read_data_from_tables(path_to_tables, 'table2_data.csv')
        self.table2_df = table2_df.set_index('layer')

    def read_data_from_tables(self, path, filename):
        fullpath = os.path.join(path, filename)
        filenameroot, file_extension = os.path.splitext(filename)
        if file_extension=='.csv':
            df = pd.read_csv(fullpath)
        elif file_extension=='.json':
                df = pd.read_json(fullpath)
        elif file_extension=='.xlsx':
                df = pd.read_excel(fullpath)
        else:
            raise NotImplementedError('Unexpected filename extension')

        return df


class Area:
    '''
    This class contains area-level data and methods
    '''
    def __init__(self, CxC, area_name='V1', requestedVFradius=.1, center_ecc=5, requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6']):
        
        self.area_name=area_name
        self.requestedVFradius=requestedVFradius
        self.center_ecc=center_ecc
        self.requested_layers=requested_layers


        # Read connection table sublayer to Table2 layer mapping. Contains assumed proportions for Ncells/sublayer
        conn_csv2table2_layer_map_df = CxC.read_data_from_tables(path_to_tables, 'connection_csv_sublayer_to_table2_layer_mapping.csv')
        self.conn_csv2table2_layer_map_df = conn_csv2table2_layer_map_df
        self.PC_apical_dendrites = CxC.read_data_from_tables(path_to_tables, 'PC_apical_dendrites.xlsx')

        # Check layer names for validity: are they mapped in the sublayer to layer mapping file
        valid_sublayers = conn_csv2table2_layer_map_df['sublayer'].tolist()
        assert set(requested_layers) - set(valid_sublayers) == set(), f'Invalid layer names, valid layer names are {valid_sublayers}'
            
        # Check that layer names match valid Table 2 layers
        assert set(conn_csv2table2_layer_map_df['layer'].tolist()) - set(CxC.table2_df.index.tolist()) == set(), f'Invalid layer to Table 2 mapping'

        if area_name=='V1':
            # Get proportion V1 of total V1
            #Requested V1 area in mm2
            V1area = self._VFradius2V1area(requestedVFradius, center_ecc)

            # Get V1 proportion for simulation so that we can get N cells from the total N cells in V1 layers in Table 2
            V1total_area = CxC.table1_df.loc['mean','V1']
            self.area_proportion = V1area / V1total_area
        else:
            self.area_proportion = 1 #Not implemented yet for other areas

        # Get layer name mapping between inhibitory cell type, excitatory cell type, Table2 and valid sublayers. You might need to update this file.
        layer_name_mapping_df = CxC.read_data_from_tables(path_to_tables, 'layer_name_mapping.xlsx')

        #Assert sublayer validity between sublayer to table2 mapping (valid_sublayers) and the current layer_name_mapping_df
        assert set(valid_sublayers) - set(layer_name_mapping_df['allowed_requested_layers'].tolist()) == set(), f'Update allowed_requested_layers in layer_name_mapping'
        self.layer_name_mapping_df = layer_name_mapping_df

    def _VFradius2V1area(self, radius, center_ecc):
        '''
        Input radius in degrees, output area in mm2
        '''
        a=1
        # Assuming M = 1/(0.077 + 0.082 × E) mm/deg Tootell 1982 Science (valid < 10 deg ecc), we get at center_ecc deg
        M = 1/(0.077 + (0.082 * center_ecc))
        # Taking 1/M = a/k + 1/k * E, we get at center_ecc deg
        k = (a + center_ecc) * M

        cx_min_mm = k * np.log(a + center_ecc - radius)
        cx_max_mm = k * np.log(a + center_ecc + radius)

        # To avoid integrating complex logarithm, we take circular patch of cx, with radius of (cx_max_mm - cx_min_mm) / 2
        radius_in_mm = (cx_max_mm - cx_min_mm) / 2
        V1area = np.pi * np.power(radius_in_mm,2)
        return V1area


class Groups:

    #TODO: AREA_OBJECT VS SELF.AREA_OBJECT STYLE

    def __init__(   self, CxC, area_object, requested_cell_types_and_proportions, cell_type_data_source, 
                    cell_type_data_folder_name, cell_type_data_file_name, monitors, bg_inputs):
        
        self.area_object = area_object
        self.requested_cell_types_and_proportions = requested_cell_types_and_proportions

        self.monitors =  monitors
        self.bg_inputs = bg_inputs

        # Unpack for init
        inhibitory_types = requested_cell_types_and_proportions['inhibitory_types']
        inhibitory_proportions = requested_cell_types_and_proportions['inhibitory_proportions']
        excitatory_types = requested_cell_types_and_proportions['excitatory_types']
        excitatory_proportions = requested_cell_types_and_proportions['excitatory_proportions']
        requested_layers = area_object.requested_layers

        # Create inhibitory_proportions_df
        # Valid EIflag 'Glutamatergic' and 'GABAergic'
        self.inhibitory_proportions_df = self.get_proportions_df(   'GABAergic',inhibitory_proportions, inhibitory_types, requested_layers, 
                                                                    cell_type_data_source, cell_type_data_folder_name, cell_type_data_file_name)
        
        self.excitatory_proportions_df = self.get_proportions_df(   'Glutamatergic',excitatory_proportions, excitatory_types, requested_layers, 
                                                                    cell_type_data_source, cell_type_data_folder_name, cell_type_data_file_name)
        
        # Choose layer mappings according to requested layers. 
        # TODO If an entry has two values separated by ";", the two values must be averaged
        layer_name_mapping_df = area_object.layer_name_mapping_df # Should be acquired from Area class instance object
        # Map requested layers to allowed neuroinformatics sublayers.
        layer_mapping_df = layer_name_mapping_df.loc[layer_name_mapping_df['allowed_requested_layers'].isin(requested_layers)]
        # Add numerical layer index to layer_mapping_df
        layer_idxs = np.arange(len(requested_layers)) + 1
        layer_mapping_df['layer_idx'] = layer_idxs.tolist()

        self.layer_mapping_df = layer_mapping_df

        # Get df with neuron groups for anatomy df, return new df to object
        self.anatomy_config_df_new = self.generate_cell_groups(CxC, area_object, requested_cell_types_and_proportions)

    def generate_cell_groups(self, CxC, area_object, requested_cell_types_and_proportions):
        '''
        generate_cell_groups function
        -for each layer in macaque V1, set the neuron groups, types and numbers
            -calculate the number of cell groups
            -generate df for holding the neuron groups
                -get starting neuron group idx and line number of insert from anatomy csv
            -map the Table2 data into cell groups
            -return df with neuron group rows for anatomy csv 
        '''
        # Unpacking for current method
        table2_df = CxC.table2_df
        anatomy_config_df = CxC.anatomy_config_df
        area_proportion = area_object.area_proportion
        requested_layers = area_object.requested_layers
        PC_apical_dendrites = area_object.PC_apical_dendrites
        inhibitory_proportions_df = self.inhibitory_proportions_df
        excitatory_proportions_df = self.excitatory_proportions_df
        monitors = self.monitors
        background_input = self.bg_inputs
        layer_mapping_df = self.layer_mapping_df 

        # Get and set column names for neuron groups
        cell_group_columns = anatomy_config_df.loc[anatomy_config_df.groupby([0]).get_group('G').index[0]-1,:].values
        existing_neuron_groups = anatomy_config_df.groupby([0]).get_group('G')
        existing_neuron_groups.columns = cell_group_columns

        # Get starting index for cell groups and row indices for anat df
        if CxC.replace_existing_cell_groups:
            start_cell_group_index = 1 # Reserve 0 for input group
            start_index = anatomy_config_df.groupby([0]).get_group('G').index[0]
        else:
            start_cell_group_index = int(existing_neuron_groups['idx'].values[-1]) + 1
            start_index = anatomy_config_df.groupby([0]).get_group('G').index[-1] + 1

        # Get the beginning of anat df before cell groups and 

        # Map cell groups to requested layers
        # Use layers as neuron subtype prefixes

        index_excitatory_s = excitatory_proportions_df.fillna(0).astype(bool).sum()
        index_inhibitory_s = inhibitory_proportions_df.fillna(0).astype(bool).sum()
        N_rows =  index_excitatory_s.sum() + index_inhibitory_s.sum()
        
        # Generate df for holding the neuron groups. Add indices here to preallocate memory
        indices = start_index + np.arange(N_rows)      
        NG_df = pd.DataFrame(columns=cell_group_columns, index=indices)
        
        # Now we go for the cell groups. First let's set the identical values
        NG_df.row_type = 'G'
        NG_df.idx = np.arange(start_cell_group_index, start_cell_group_index + N_rows)
        NG_df.net_center = '--'
        NG_df.monitors = monitors

        # Add one layer at a time starting from L1
        current_group_index = start_index

        for layer in requested_layers:
            current_layer_proportion_inhibitory = table2_df.loc[layer]['percent_inhibitory'] / 100
            current_layer_N_neurons_area = table2_df.loc[layer]['n_neurons_10e6'] * 1e6
            current_layer_N_excitatory_neurons = (1 - current_layer_proportion_inhibitory) * area_proportion * current_layer_N_neurons_area
            current_layer_N_inhibitory_neurons = current_layer_proportion_inhibitory * area_proportion * current_layer_N_neurons_area
            layer_idx = layer_mapping_df.loc[layer_mapping_df['allowed_requested_layers']==layer]['layer_idx'].values

            # Add excitatory cell groups
            for current_group in excitatory_proportions_df.index.values:
                # Jump to next group if proportion is 0
                if excitatory_proportions_df.loc[current_group][layer]==0:
                    continue
                # go through all columns which require individual values
                NG_df.loc[current_group_index,'number_of_neurons'] =  np.round(excitatory_proportions_df.loc[current_group][layer] * current_layer_N_excitatory_neurons)
                NG_df.loc[current_group_index,'neuron_type'] = current_group
                NG_df.loc[current_group_index,'neuron_subtype'] = layer + '_' + current_group
                # For PC, map apical dendrite extent from table to layer index
                if current_group == 'PC':
                    # TODO: NOT IMPLEMENTED YET
                    pass
                    # NG_df.loc[current_group_index,layer_idx] =
                else: 
                    NG_df.loc[current_group_index,'layer_idx'] =  layer_idx

                NG_df.loc[current_group_index,'n_background_inputs'] = background_input['n_background_inputs_for_excitatory_neurons']
                NG_df.loc[current_group_index,'n_background_inhibition'] = background_input['n_background_inhibition_for_excitatory_neurons']
                current_group_index += 1

            # Add inhibitory cell groups
            for current_group in inhibitory_proportions_df.index.values:
                # Jump to next group if proportion is 0
                if inhibitory_proportions_df.loc[current_group][layer]==0:
                    continue
                # go through all columns which require individual values
                NG_df.loc[current_group_index,'number_of_neurons'] =  np.round(inhibitory_proportions_df.loc[current_group][layer] * current_layer_N_inhibitory_neurons)
                NG_df.loc[current_group_index,'neuron_type'] = current_group
                NG_df.loc[current_group_index,'neuron_subtype'] = layer + '_' + current_group
                
                NG_df.loc[current_group_index,'layer_idx'] =  layer_idx

                NG_df.loc[current_group_index,'n_background_inputs'] = background_input['n_background_inputs_for_inhibitory_neurons']
                NG_df.loc[current_group_index,'n_background_inhibition'] = background_input['n_background_inhibition_for_inhibitory_neurons']
                current_group_index += 1

        assert current_group_index == indices[-1] + 1, 'oh-ou'

        # Add to anatomy df

        # Change column names
        NG_df.columns = anatomy_config_df.columns

        # Get end of cell groups index for slicing original anat df
        end_index = anatomy_config_df.groupby([0]).get_group('G').index[-1] + 1
        anatomy_config_df_beginning = anatomy_config_df.iloc[:start_index,:] 
        anatomy_config_df_end = anatomy_config_df.iloc[end_index:,:] 
        anatomy_config_df_new = pd.concat([anatomy_config_df_beginning, NG_df, anatomy_config_df_end], ignore_index=True)

        return anatomy_config_df_new

    def drop_unused_cell_types(self, proportions_df, types):
            unused_types_set = set(proportions_df.index) - set(types)
            proportions_df_clean = proportions_df.drop(unused_types_set, axis=0)
            # rescale to sum 1 after dropping the unused types
            proportions_df = proportions_df_clean / proportions_df_clean.sum(axis=0)
            return proportions_df

    def get_markram_cell_type_proportions(self, cell_type_df, EIflag):
        # This diverges according to EI flag because excitatory is not implemented currently
        
        if EIflag=='Glutamatergic':
            # This can be fixed later. Now it gives proportions of cell types as 1/N types
            excitatory_proportions_df = self.get_empty_cell_type_proportions(self.area_object.requested_layers, self.requested_cell_types_and_proportions['excitatory_types'])
            return excitatory_proportions_df

        elif EIflag=='GABAergic':
            # Get cell type proportions from Markram rat somatosensory cx data
            n_neurons_per_type = cell_type_df.loc['No. of neurons per morphological types',:]
            type_mapping = {
                'SST' : ['MC'],
                'PVALB' : ['NBC','LBC'],
                'VIP' : ['SBC','BP','DBC']}
            layers_for_cell_count = cell_type_df.columns.values
            type_count = {}
            inhibitory_proportions_df = pd.DataFrame(index=type_mapping.keys())

            # calculate sum of inhibitory SST, PV and VIP cells
            for layer_for_count in layers_for_cell_count:
                cell_count_dict = n_neurons_per_type[layer_for_count]
                layer_string_to_add = layer_for_count + '_'
                for current_type in type_mapping.keys():
                    n_neurons_for_current_type = 0
                    subtypes_list = type_mapping[current_type]
                    for current_subtype in subtypes_list:
                        try:
                            n_neurons_for_current_subtype = cell_count_dict[layer_string_to_add + current_subtype]
                            n_neurons_for_current_type += n_neurons_for_current_subtype
                        except KeyError:
                            continue
                    type_count[current_type] = n_neurons_for_current_type

                # get proportions of the three types in different layers
                type_count_array = np.fromiter(type_count.values(), dtype=float)
                proportions = type_count_array/sum(type_count_array)
                inhibitory_proportions_df[layer_for_count] = pd.Series(proportions, index=type_mapping.keys())
            return inhibitory_proportions_df

    def get_allen_cell_type_proportions(self, cell_type_allen_df_raw, EIflag):
        # Valid EIflag 'Glutamatergic' and 'GABAergic'

        # Select appropriate classifying parameters from df
        cell_type_allen_relevant_columns_df = cell_type_allen_df_raw[['region_label', 'class_label', 'subclass_label', 'cortical_layer_label']]
        cell_type_allen_V1_df = cell_type_allen_relevant_columns_df.loc[cell_type_allen_relevant_columns_df['region_label'] == 'V1C']
        cell_type_allen_df = cell_type_allen_V1_df.drop(columns=['region_label'])
        
        grouped = cell_type_allen_df.groupby(['class_label', 'cortical_layer_label'])
        allen_layers = cell_type_allen_df['cortical_layer_label'].unique()

        # Prepare for getting unique cell types from Allen data
        cell_type_allen_df_for_unique_types = cell_type_allen_df.drop('cortical_layer_label', axis=1)
        grouped_for_unique_types = cell_type_allen_df_for_unique_types.groupby(['class_label'])

        # Neuron proportions
        # Get unique types
        types=grouped_for_unique_types.get_group(EIflag)['subclass_label'].unique()
        proportions_df = pd.DataFrame(index=types)

        for layer_for_count in allen_layers:
            neurons_by_layer = grouped.get_group((EIflag,layer_for_count))
            proportions = neurons_by_layer['subclass_label'].value_counts(normalize=True)
            proportions_df[layer_for_count] = proportions

        return proportions_df

    def get_empty_cell_type_proportions(self, requested_layers, requested_types):
        '''
        Create cell type proportions from N cell types when they are not defined by user or pre-existing data
        '''
        N_types = len(requested_types)
        proportions = 1/N_types
        proportions_df = pd.DataFrame(index=requested_types, columns=requested_layers).fillna(proportions)

        return proportions_df

    def get_proportions_df(self, EIflag, proportions, types, requested_layers, cell_type_data_source, cell_type_data_folder_name, cell_type_data_file_name):


        if proportions: # If manually given, executes this end continues
            proportions_df =  pd.DataFrame(data=proportions, index=types)

        elif cell_type_data_source: # If given
            fullpath = os.path.join(path_to_tables, cell_type_data_folder_name)
            cell_type_df = CxC.read_data_from_tables(fullpath, cell_type_data_file_name)

            if cell_type_data_source=='HBP' or cell_type_data_source=='Markram':
                proportions_df = self.get_markram_cell_type_proportions(cell_type_df, EIflag)

            elif cell_type_data_source=='Allen':
                # Valid EIflag 'Glutamatergic' and 'GABAergic'
                proportions_df = self.get_allen_cell_type_proportions(cell_type_df, EIflag)

            proportions_df = self.drop_unused_cell_types(proportions_df, types)

        # If no proportions is defined manually and no data source is given, fill proportions with 1/N cell types
        elif not cell_type_data_source: 
            proportions_df = self.get_empty_cell_type_proportions(requested_layers, types)

        return proportions_df


class Connections:

    def __init__(self, CxC, NG_new):
        '''
        Generate synapses object, which includes the new anatomy df with connections
        '''
        
        # Set parameters to connections object
        self.replace_existing_cell_groups = CxC.replace_existing_cell_groups
        
        # Read data from files.
        # Read ni csv into dataframe
        self.exc_df = CxC.read_data_from_tables(path_to_ni_csv, 'connections_local_excitatory.csv')
        self.inh_df = CxC.read_data_from_tables(path_to_ni_csv, 'connections_local_inhibitory.csv')

        self.anatomy_config_df_new_connections = NG_new.anatomy_config_df_new

        self.anatomy_config_df_new_connections_new_synapses = self.generate_synapses()

        pdb.set_trace()


    def generate_synapses(self):
        '''
        generate_synapses function
        -for each origin and target group, set receptor,pre_syn_idx,post_syn_idx by layer and compartment,syn_type,p,n,
        -return df with receptor,pre_syn_idx,post_syn_idx,syn_type,p,n, 

        '''

        # Unpack for this method
        replace_existing_cell_groups = self.replace_existing_cell_groups
        exc_df = self.exc_df 
        inh_df = self.inh_df 
        anatomy_config_df_new_connections = self.anatomy_config_df_new_connections

        pdb.set_trace()

        pass

def show_cortex():
    pass
'''
show_cortex
Visualize the neurons (and connections?) in 3D
visimpl?
'''


    

'''
layers = [L1, L23, L4A,L4B, L4CΑ, L4CB,L5,L6]; We have the N cells for these layers in Table 2
'''
'''
ni2anat function
A filter for turning neuroinfo csv data into connections in anatomical csv file
-read ni csv into dataframe
-input layers,  
    
The sublayer to layer mapping file contains eg sublayers 4Cm, L6A etc. Cell proportions in sublayers are simply
1/(N sub layers). The blob/interblob relation is estimated to be 0.4/0.6, according to surface areas at 5 deg ecc (Livingstone 1984 JNeurosci)

-check their names for validity from the sublayer to layer mapping file

-read in Table 1 and Table 2 data
-map the N layers to valid Table 2 layers, 

-call generate_cell_groups to get df with neuron group names, numbers & positions by layer and type
-call cell_group_row for anat csv row output for cell groups

-map each row in the ni csv file into layer of origin and termination, set strength
-set origin groups according to subgroups
-set target groups according to subgroups
-call generate_synapses function
-call synapse_row

-read in anat csv to start with
-write cell groups and synapses into the existing anat csv file

-call show_cortex

'''

if __name__ == "__main__":
    '''
    Start of user input
    Copy and comment/uncomment examples/your own versions by need. If python gives exception, look first your own syntax below.
    '''

    area_name='V1' # Don't change this.
    requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6'] # You should be able start from L4CA or B alone for testing
    requestedVFradius=.1 # Increasing this rapidly makes much more cells and increases the computational cost
    center_ecc=5 # Don't change this. This might be later replaced by 2D coordinates

    '''
    The proportion of inhibitory and excitatory neurons in distinct V1 layers will come from our review Table2.
    Below, you provide inhibitory and excitatory cell types and proportions when you have more than one type for each.
    Cell types are mandatory, and their names must match both Allen/HBP data table if these are used and physiology file cell type names
    Inhibitory and excitatory cell type proportions layerwise come either from  Allen/HBP data, 
    or you can define them by hand below. If left empty, 1/N cell types will be used for each layer.

    Later in an advanced version of the system, you can use the considerations below:
    Inhibitory neurons are divided into three groups in each layer, except in L1. In L1 we have limited data. The inh neurons will
    for sure be structurally simple. Whether we need a separate cell type for L1 is an open question. We could use the LAMP5 molecular marker as L1I 
    neuron, but then we need to assing its proportion in other layers with no idea of its structure. See Tasic_2018_Nature
    One limitation is that we are using Allen molecular types to get quantities of our limited set of structural types.
    inhibitory_types = ['SST', 'VIP', 'PVALB'] # This is a good guess, and take these proportions from Allen data

    For excitatory neurons we do not have HBP/Allen data to start with because only structural types currently implementd in CxSystem are 
    pyramidal cells (PC) with apical dendrites and spiny stellate (SS) cells which are point-like.
    Each PC group in layers L2 to L6 are divided into two. The first has apical dendrite extending to L1 and the second to L23 (L4-5, or L4C (L6)
    see Lund_1991_JCompNeurol, Callaway_1996_VisNeurosci, Nassi_2007_Neuron, Nhan_2012_JCompNeurol, Wiser_1996_Jneurosci, Briggs_2016_Neuron

    '''
    
    # Here are some examples
    # inhibitory_types = ['SST', 'VIP', 'PVALB']
    inhibitory_types = ['PVALB']
    inhibitory_proportions={} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
    # inhibitory_proportions = {  
    # 'L1': [1, 0, 0], 
    # 'L23': [0, .5, .5], 
    # 'L4A': [.33, .33, .33], 
    # 'L4B': [.33, .33, .33], 
    # 'L4C': [1, 0, 0],
    # 'L5': [0, .5, .5], 
    # 'L6': [0, .5, .5]}

    # Excitatory proportions are given by hand here. 
    # The list length for each layer must match the N types and should sum to approximately 1.
    # The layers must match the requested layers
    # The SS type in L1 is just a pointlike excitatory neuron
    excitatory_types = ['SS', 'PC1', 'PC2']
    # excitatory_types = ['SS'] # IT is one of th Allen excitatory types. SS = spiny stellate and this should exist in all physiological files
    # excitatory_proportions = {} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
    excitatory_proportions = {  
    'L1': [1, 0, 0], 
    'L23': [0, .5, .5], 
    'L4A': [.33, .33, .33], 
    'L4B': [.33, .33, .33], 
    'L4CA': [1, 0, 0],
    'L4CB': [1, 0, 0],
    'L5': [0, .5, .5], 
    'L6': [0, .5, .5]}

    # inhibitory_types = ['example1','example2','example3'] # Name example is used for construction. For simulation, these need to match the physiology cell type names
    # inhibitory_proportions={}

    # excitatory_types = ['example4'] 
    # excitatory_proportions = {}

    # Read in anat csv to start with. Check for config_files folder for valid inputs.
    # If replace_existing_cell_groups flag is False, your groups will be added to current groups. This allows building system stepwise
    replace_existing_cell_groups = True # 
    anatomy_config_file_name = 'pytest_anatomy_config.csv'
    physiology_config_file_name = 'pytest_physiology_config.csv' # anatomy and physiology filenames sometimes diverge

    
    # Activate correct selection. These data provide proportions of distinct cell groups in each layer.
    # Markram data, rat S1
    # cell_type_data_source = 'HBP'
    
    # # Allen data, human V1
    # cell_type_data_source = 'Allen'

    # No data source. Use this eg for single excitatory and single inhibitory types. 
    cell_type_data_source = ''

    request_monitors = '[Sp]'
    n_background_inputs_for_excitatory_neurons = 630
    n_background_inhibition_for_excitatory_neurons = 290    
    n_background_inputs_for_inhibitory_neurons = 500
    n_background_inhibition_for_inhibitory_neurons = 180    

    '''
    End of user input
    '''
    # TODO: ASSERT excitatory_proportions IS EITHER EMPTY, OR ROWS MATCHES N NEURON TYPES AND COLUMN NAMES MATCHES THE REQUESTED LAYERS

    if cell_type_data_source == 'HBP':
        cell_type_data_folder_name='hbp_data'; cell_type_data_file_name='layer_download.json'
    elif cell_type_data_source == 'Allen':
        cell_type_data_folder_name='allen_data'; cell_type_data_file_name='sample_annotations.csv'
    elif cell_type_data_source == '':
        cell_type_data_folder_name=''; cell_type_data_file_name=''

    # Packing of variables for brevity
    requested_cell_types_and_proportions = {
        'inhibitory_types':inhibitory_types, 
        'inhibitory_proportions':inhibitory_proportions, 
        'excitatory_types':excitatory_types,
        'excitatory_proportions':excitatory_proportions}
    requested_background_input = {
        'n_background_inputs_for_excitatory_neurons':n_background_inputs_for_excitatory_neurons,
        'n_background_inhibition_for_excitatory_neurons':n_background_inhibition_for_excitatory_neurons,
        'n_background_inputs_for_inhibitory_neurons':n_background_inputs_for_inhibitory_neurons,
        'n_background_inhibition_for_inhibitory_neurons':n_background_inhibition_for_inhibitory_neurons}

    anatomy_config_df = read_config_file(os.path.join(path_to_config_files,anatomy_config_file_name))
    physiology_config_df = read_config_file(os.path.join(path_to_config_files,physiology_config_file_name))    

    CxC = Config(replace_existing_cell_groups=replace_existing_cell_groups, anatomy_config_df=anatomy_config_df, physiology_config_df=physiology_config_df)
    V1 = Area(CxC, area_name=area_name, requestedVFradius=requestedVFradius, center_ecc=center_ecc, requested_layers=requested_layers)

    # Add anatomy and physiology config files to start with
    suffix_for_new_file = '_cxc'

    NG_new = Groups(  CxC, V1, requested_cell_types_and_proportions, cell_type_data_source, cell_type_data_folder_name, 
                cell_type_data_file_name, request_monitors,requested_background_input)

    NG_new.anatomy_config_df_new.to_csv(os.path.join(path_to_config_files,anatomy_config_file_name[:-4] + suffix_for_new_file + '.csv'), header=False, index=False)

    Conn_new = Connections(CxC, NG_new)

