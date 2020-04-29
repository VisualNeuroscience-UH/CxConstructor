import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd
import pdb


path_to_tables = r'C:\Users\Simo\Laskenta\Models\MacV1Buildup\tables'
path_to_ni_csv = r'C:\Users\Simo\Laskenta\Models\MacV1Buildup\ni_csv_copy'

# Module level utilities
def read_data_from_tables(path, filename):
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
    def __init__(self, area_name='V1', requestedVFradius=.1, center_ecc=5, requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6']):
        
        self.area_name=area_name
        self.requestedVFradius=requestedVFradius
        self.center_ecc=center_ecc
        self.requested_layers=requested_layers


        # Read connection table sublayer to Table2 layer mapping. Contains assumed proportions for Ncells/sublayer
        sublayers_df = read_data_from_tables(path_to_tables, 'connection_csv_sublayer_to_table2_layer_mapping.csv')
        self.sublayers_df = sublayers_df

        # Check layer names for validity: are they mapped in the sublayer to layer mapping file
        valid_sublayers = sublayers_df['sublayer'].tolist()
        assert set(requested_layers) - set(valid_sublayers) == set(), f'Invalid layer names, valid layer names are {valid_sublayers}'
        
        # Read data from Table 1 (surface area in mm2 for V1, V2 and V5) and Table 2 (N neurons and synapses for V1 in distinct layers) 
        table1_df = read_data_from_tables(path_to_tables, 'table1_data.csv')
        self.table1_df = table1_df.set_index('stat') 
        table2_df = read_data_from_tables(path_to_tables, 'table2_data.csv')
    
        # Check that layer names match valid Table 2 layers
        assert set(sublayers_df['layer'].tolist()) - set(table2_df['layer'].tolist()) == set(), f'Invalid layer to Table 2 mapping'
        self.table2_df = table2_df

        if area_name=='V1':
            # Get proportion V1 of total V1
            #Requested V1 area in mm2
            V1area = self._VFradius2V1area(requestedVFradius, center_ecc)

            # Get V1 proportion for simulation so that we can get N cells from the total N cells in V1 layers in Table 2
            V1total_area = self.table1_df.loc['mean','V1']
            self.area_proportion = V1area / V1total_area
        else:
            self.area_proportion = 1 #Not implemented yet for other areas

        # Get layer name mapping between inhibitory cell type, excitatory cell type, Table2 and valid sublayers. You might need to update this file.
        layer_name_mapping_df = read_data_from_tables(path_to_tables, 'layer_name_mapping.xlsx')

        #Assert sublayer validity between sublayer to table2 mapping (valid_sublayers) and the current layer_name_mapping_df
        assert set(valid_sublayers) - set(layer_name_mapping_df['requested_layers'].tolist()) == set(), f'Update requested_layers in layer_name_mapping'
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


class Groups():

    def __init__(self, area_object, requested_cell_types_and_proportions, cell_type_data_folder_name, cell_type_data_file_name, cell_type_data_source='Allen'):
        
        fullpath = os.path.join(path_to_tables, cell_type_data_folder_name)
        cell_type_df = read_data_from_tables(fullpath, cell_type_data_file_name)

        if cell_type_data_source=='HBP' or cell_type_data_source=='Markram':
            # Get inhibitory cell type proportions from Markram rat somatosensory data
            self.inhibitory_proportions_df = get_markram_cell_type_proportions(cell_type_df)

        elif cell_type_data_source=='Allen':
            # Get inhibitory cell type proportions from Allen institute data for human V1
            # Select appropriate classifying parameters from df
            cell_type_allen_relevant_columns_df = cell_type_df[['region_label', 'class_label', 'subclass_label', 'cortical_layer_label']]
            # Get unique labels for comparison
            # valid_regions=set(cell_type_allen_relevant_columns_df['region_label']) # In case u want to compare regions later
            # valid_layers=set(cell_type_allen_relevant_columns_df['cortical_layer_label']) # In case u want to compare regions later
            # valid_class_labels=set(cell_type_allen_relevant_columns_df['class_label']) # In case u want to compare regions later
            # valid_subclass_labels=set(cell_type_allen_relevant_columns_df['subclass_label']) # In case u want to compare regions later
            # Select visual cortex according to region label
            cell_type_allen_V1_df = cell_type_allen_relevant_columns_df.loc[cell_type_allen_relevant_columns_df['region_label'] == 'V1C']
            cell_type_allen_V1clean_df = cell_type_allen_V1_df.drop(columns=['region_label'])
        
            # Excitatory proportions are not currently used
            self.inhibitory_proportions_df, self.excitatory_proportions_df = get_allen_cell_type_proportions(cell_type_allen_V1clean_df)

        # Get df with neuron group names, numbers & positions by layer and type
        NG_df = self.generate_cell_groups(area_object, requested_cell_types_and_proportions)

    def generate_cell_groups(self, area_object, requested_cell_types_and_proportions):
        '''
        generate_cell_groups function
        -for each layer in macaque V1, set the neuron groups, types and numbers
            -unit N / mm2 / layer
            -input 
                -area_proportion: proportion of V1 area for simulation, scales N neurons from Table 2
                -requested_layers: requested layers for constructing anatomy csv for simulations
                -table2_df: Vanni et al 2020 CerCx Table2 data, ie layerwise  N neurons and proportion of inhibitory cells
                -inhibitory_proportions_allen_df: Human data (allen) or rodent data (HBP/Markram) for layerwise proportions of neuron types
            -set inhibitory types
            -calculate the number of cell groups
            -generate df for holding the neuron groups
            row_type,idx,number_of_neurons,neuron_type,neuron_subtype,layer_idx,net_center,monitors,n_background_inputs,n_background_inhibition
            -map the Table2 data into cell groups
            -calculate density for each group
            -set positions by layer and type
            -return df with neuron group names, numbers & positions by layer and type 
            TÄHÄN JÄIT: LAITA INH TYPES, EXC TYPES JA EXC PROPORTIONS JÄRKEVÄMPÄÄN PAIKKAAN, JONNEKIN ESILLE
        '''
        # Unpacking for current method
        area_proportion = area_object.area_proportion
        requested_layers = area_object.requested_layers
        table2_df = area_object.table2_df
        inhibitory_proportions_df = self.inhibitory_proportions_df
        excitatory_proportions_df = self.excitatory_proportions_df
        inhibitory_types = requested_cell_types_and_proportions['inhibitory_types']
        inhibitory_proportions = requested_cell_types_and_proportions['inhibitory_proportions']
        excitatory_types = requested_cell_types_and_proportions['excitatory_types']
        excitatory_proportions = requested_cell_types_and_proportions['excitatory_proportions']

        # Calculate the number of cell groups

        # Set inhibitory types
        # inhibitory neurons are divided into three groups in each layer, except in L1. In L1 we have limited data. The inh neurons will
        # for sure be structurally simple. Whether we need a separate cell type is an open question. We could use the LAMP5 molecular marker as L1I 
        # neuron, but then we need to assing its proportion in other layers with no idea of its structure. See Tasic_2018_Nature
        # One limitation is that we are using Allen molecular types to get quantities of our limited set of structural types.
        # inhibitory_types = ['SST', 'VIP', 'PVALB']

        if inhibitory_proportions:
            inhibitory_proportions_df =  pd.DataFrame(data=inhibitory_proportions, index=inhibitory_types)
        else:
            unused_types_set = set(inhibitory_proportions_df.index) - set(inhibitory_types)
            inhibitory_proportions_df_clean = inhibitory_proportions_df.drop(unused_types_set, axis=0)
            # rescale to sum 1 after dropping the unused types
            inhibitory_proportions_df = inhibitory_proportions_df_clean / inhibitory_proportions_df_clean.sum(axis=0)

        # # Set excitatory types
        # # Each PC group in layers L2 to L6 are divided into two. The first has apical dendrite extending to L1 and the second to L23 (L4-5, or L4C (L6)
        # # see Lund_1991_JCompNeurol, Callaway_1996_VisNeurosci, Nassi_2007_Neuron, Nhan_2012_JCompNeurol, Wiser_1996_Jneurosci, Briggs_2016_Neuron
        # # SS groups are in all L4 sublayers
        # excitatory_types = ['SS', 'PC1', 'PC2']
        # # Construct df manually, excitatory_proportions should sum to about 1
        # excitatory_proportions = {  
        #     'L1': [0, 1, 0], 
        #     'L23': [0, .5, .5], 
        #     'L4A': [.33, .33, .33], 
        #     'L4B': [.33, .33, .33], 
        #     'L4C': [1, 0, 0],
        #     'L5': [0, .5, .5], 
        #     'L6': [0, .5, .5]}
        if excitatory_proportions:
            excitatory_proportions_df = pd.DataFrame(data=excitatory_proportions, index=excitatory_types)
        else:
            unused_types_set = set(excitatory_proportions_df.index) - set(excitatory_types)
            excitatory_proportions_df_clean = excitatory_proportions_df.drop(unused_types_set, axis=0)
            # rescale to sum 1 after dropping the unused types
            excitatory_proportions_df = excitatory_proportions_df_clean / excitatory_proportions_df_clean.sum(axis=0)

        # Choose layer mappings according to requested layers. If an entry has two values separated by ";", the two values must be averaged
        layer_name_mapping_df = area_object.layer_name_mapping_df # Should be inherited from Area class
        layer_mapping_df = layer_name_mapping_df.loc[layer_name_mapping_df['requested_layers'].isin(requested_layers)]

        # TÄHÄN JÄIT: ASSERTIOT E JA I DFLLE, TESTAA KAIKKI NELJÄ IF KONDITIOTA
        # Generate df for holding the neuron groups
        columns = ['row_type','idx','number_of_neurons','neuron_type','neuron_subtype','layer_idx','net_center','monitors','n_background_inputs','n_background_inhibition']
        NG_df = pd.DataFrame(columns=columns)

        # Now we go for the cell groups

        pdb.set_trace()

    def generate_cell_rows(self):
            # Generate df for holding the neuron groups
        columns = ['row_type','idx','number_of_neurons','neuron_type','neuron_subtype','layer_idx','net_center','monitors','n_background_inputs','n_background_inhibition']
        cell_rows_df = pd.DataFrame(columns=columns)
        #


class Connections:

    def __init__(self):
        
        # Read data from files.
        # Read ni csv into dataframe
        self.exc_df = read_data_from_tables(path_to_ni_csv, 'connections_local_excitatory.csv')
        self.inh_df = read_data_from_tables(path_to_ni_csv, 'connections_local_inhibitory.csv')



'''
How do we assign connections to distinct cell groups -- Peter's law? 
-- start from Peters' rule. Later, includes White's exceptions (see Braitenberg & Schuz 1998, pp 100-101)
    -create cell numbers and type proportions by layer -- lue
        -excitatory: SS and PC in L4, PC elsewhere
            -excitatory subgroups?
        -inhibitory: PV, SST and VIP neurons
            -proportions in distinct layers, data for macaques? If not then rodents -- eg Markram

Make up naming and coding system for cell groups in distinct layers according to CxSystem framework
    -preparation for flexibility
 Table 1: Cortical surface areas (mm2) from anatomical studies
 Table 2: Total number of neurons, synapses/neuron, and the proportion of inhibitory interneurons in each cortical layer of area V1
 Table 4: Summary of horizontal connectivity
'''


def get_markram_cell_type_proportions(cell_type_df):
    # Get cell type proportions from Markram rat somatosensory data
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

def get_allen_cell_type_proportions(cell_type_allen_df):
    grouped = cell_type_allen_df.groupby(['class_label', 'cortical_layer_label'])
    allen_layers = cell_type_allen_df['cortical_layer_label'].unique()
    # n_allen_layers = cell_type_allen_df['cortical_layer_label'].nunique()

    # Prepare for getting unique cell types from Allen data
    cell_type_allen_df_for_unique_types = cell_type_allen_df.drop('cortical_layer_label', axis=1)
    grouped_for_unique_types = cell_type_allen_df_for_unique_types.groupby(['class_label'])


    # # Inhibitory proportions
    # types=['SST','PVALB','VIP']
    # #  Empty df for collecting inhibitory neuron proportions by layer
    # inhibitory_proportions_df = pd.DataFrame(index=types)

    # for layer_for_count in allen_layers:
    #     GABAneurons_by_layer_df = grouped.get_group(('GABAergic',layer_for_count))
    #     mask = GABAneurons_by_layer_df.subclass_label.isin(types)
    #     GABAneurons_by_layer_df_clean = GABAneurons_by_layer_df[mask]
    #     proportions = GABAneurons_by_layer_df_clean['subclass_label'].value_counts(normalize=True)
    #     inhibitory_proportions_df[layer_for_count] = proportions

    # Excitatory proportions
    # Get unique types
    types_exc=grouped_for_unique_types.get_group('Glutamatergic')['subclass_label'].unique()
    types_inh=grouped_for_unique_types.get_group('GABAergic')['subclass_label'].unique()
    excitatory_proportions_df = pd.DataFrame(index=types_exc)
    inhibitory_proportions_df = pd.DataFrame(index=types_inh)

    for layer_for_count in allen_layers:
        GLUneurons_by_layer = grouped.get_group(('Glutamatergic',layer_for_count))
        GABAneurons_by_layer = grouped.get_group(('GABAergic',layer_for_count))
        proportions_exc = GLUneurons_by_layer['subclass_label'].value_counts(normalize=True)
        proportions_inh = GABAneurons_by_layer['subclass_label'].value_counts(normalize=True)
        excitatory_proportions_df[layer_for_count] = proportions_exc
        inhibitory_proportions_df[layer_for_count] = proportions_inh

    return inhibitory_proportions_df, excitatory_proportions_df
 
def cell_group_row():
    '''
    cell_group_row function
    We want to keep this separate method, because the syntax of the output csv file may change
    -return row_type,idx,number_of_neurons,neuron_type,neuron_subtype,layer_idx,net_center,monitors,n_background_inputs,n_background_inhibition
    '''
    pass

def generate_synapses():
    '''
    generate_synapses function
    -for each origin and target group, set receptor,pre_syn_idx,post_syn_idx by layer and compartment,syn_type,p,n,
    -return df with receptor,pre_syn_idx,post_syn_idx,syn_type,p,n, 
    '''
    pass

def synapse_row():
    pass
'''
synapse_row function
We want to keep this separate method, because the syntax of the output csv file may change
-return df with: row_type: "S" ,receptor,pre_syn_idx,post_syn_idx,syn_type,p,n,monitors,load_connection,save_connection,custom_weight,spatial_decay
'''

def show_cortex():
    pass
'''
show_cortex
Visualize the neurons (and connections?) in 3D
visimpl?
'''

# def ni2anat(requestedVFradius=.1, center_ecc=5,layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6']):



    

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
    
    area_name='V1'
    requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6']
    requestedVFradius=.1
    center_ecc=5

    # Inhibitory proportions come from Allen data below
    # Their types must match the cell_type_data_file_name data
    inhibitory_types = ['SST', 'VIP', 'PVALB']
    inhibitory_proportions={}+-

    # Excitatory types are given by hand here
    excitatory_types = ['SS', 'PC1', 'PC2']
    excitatory_proportions = {  
    'L1': [0, 1, 0], 
    'L23': [0, .5, .5], 
    'L4A': [.33, .33, .33], 
    'L4B': [.33, .33, .33], 
    'L4C': [1, 0, 0],
    'L5': [0, .5, .5], 
    'L6': [0, .5, .5]}

    # # Inhibitory proportions come from Allen data below
    # # Their types must match the cell_type_data_file_name data
    # inhibitory_types = ['PVALB']
    # inhibitory_proportions={}

    # # Excitatory types are given by hand here
    # excitatory_types = ['IT'] # This IT is an excitatory cell type in Allen data, referes to intratelencephalic
    # excitatory_proportions = {}

    # Let's pack them for brevity
    requested_cell_types_and_proportions = {
        'inhibitory_types':inhibitory_types, 
        'inhibitory_proportions':inhibitory_proportions, 
        'excitatory_types':excitatory_types,
        'excitatory_proportions':excitatory_proportions}

    V1=Area(area_name=area_name, requestedVFradius=requestedVFradius, center_ecc=center_ecc, requested_layers=requested_layers)
    
    # # Markram data, rat S1
    # cell_type_data_folder_name='hbp_data'; cell_type_data_file_name='layer_download.json'
    
    # Allen data, human V1
    cell_type_data_folder_name='allen_data'; cell_type_data_file_name='sample_annotations.csv'

    xx=Groups(V1, requested_cell_types_and_proportions, cell_type_data_folder_name, cell_type_data_file_name, cell_type_data_source='Allen')