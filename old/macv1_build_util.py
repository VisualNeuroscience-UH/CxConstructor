import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd
import pdb


path_to_tables = r'C:\Users\Simo\Laskenta\Models\MacV1Buildup\tables'
path_to_ni_csv = r'C:\Users\Simo\Laskenta\Models\MacV1Buildup\ni_csv_copy'


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
def generate_cell_groups(v1_proportion, requested_layers, table2_df, inhibitory_proportions_allen_df, layer_name_mapping_df):
    '''
    generate_cell_groups function
    -for each layer in macaque V1, set the neuron groups, types and numbers
        -unit N / mm2 / layer
        -input 
            -v1_proportion: proportion of V1 area for simulation, scales N neurons from Table 2
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
    '''
    # Calculate the number of cell groups

    # Set inhibitory types
    # inhibitory neurons are divided into three groups in each layer, except in L1. In L1 we have limited data. The inh neurons will
    # for sure be structurally simple. Whether we need a separate cell type is an open question. We could use the LAMP5 molecular marker as L1I 
    # neuron, but then we need to assing its proportion in other layers with no idea of its structure. See Tasic_2018_Nature
    # One limitation is that we are using Allen molecular types to get quantities of our limited set of structural types.
    inhibitory_types = ['SST', 'VIP', 'PVALB']
    unused_types_set = set(inhibitory_proportions_allen_df.index) - set(inhibitory_types)
    inhibitory_proportions_allen_df_clean = inhibitory_proportions_allen_df.drop(unused_types_set, axis=0)
    # rescale to sum 1 after dropping the unused types
    inhibitory_proportions_df = inhibitory_proportions_allen_df_clean / inhibitory_proportions_allen_df_clean.sum(axis=0)

    # Set excitatory types
    # Each PC group in layers L2 to L6 are divided into two. The first has apical dendrite extending to L1 and the second to L23 (L4-5, or L4C (L6)
    # see Lund_1991_JCompNeurol, Callaway_1996_VisNeurosci, Nassi_2007_Neuron, Nhan_2012_JCompNeurol, Wiser_1996_Jneurosci, Briggs_2016_Neuron
    # SS groups are in all L4 sublayers
    excitatory_types = ['SS', 'PC1', 'PC2']
    # Construct df manually, excitatory_proportions should sum to about 1
    excitatory_proportions = {  
        'L1': [0, 1, 0], 
        'L23': [0, .5, .5], 
        'L4A': [.33, .33, .33], 
        'L4B': [.33, .33, .33], 
        'L4C': [1, 0, 0],
        'L5': [0, .5, .5], 
        'L6': [0, .5, .5]}
    excitatory_proportions_df = pd.DataFrame(data=excitatory_proportions, index=excitatory_types)

    # Choose layer mappings according to requested layers. If an entry has two values separated by ";", the two values must be averaged
    layer_mapping_df = layer_name_mapping_df.loc[layer_name_mapping_df['requested_layers'].isin(requested_layers)]

    # Generate df for holding the neuron groups
    columns = ['row_type','idx','number_of_neurons','neuron_type','neuron_subtype','layer_idx','net_center','monitors','n_background_inputs','n_background_inhibition']
    NG_df = pd.DataFrame(columns=columns)

    # Now we go for the cell groups

    pdb.set_trace()
    pass

def VFradius2V1area(radius, center_ecc):
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

def ni2anat(requestedVFradius=.1, center_ecc=5,layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6']):
    # This function reads all data from files.
    # Read ni csv into dataframe
    excitatory_connections_filename = 'connections_local_excitatory.csv'
    inhibitory_connections_filename = 'connections_local_inhibitory.csv'
    exc_fullpath = os.path.join(path_to_ni_csv, excitatory_connections_filename)
    inh_fullpath = os.path.join(path_to_ni_csv, inhibitory_connections_filename)
    exc_df = pd.read_csv(exc_fullpath)
    inh_df = pd.read_csv(inh_fullpath)

    # Read connection table sublayer to Table2 layer mapping
    sublayer_to_layer_mapping_filename = 'connection_csv_sublayer_to_table2_layer_mapping.csv'
    sublayers_fullpath = os.path.join(path_to_tables, sublayer_to_layer_mapping_filename)
    sublayers_df = pd.read_csv(sublayers_fullpath)

    # Check layer names for validity: are they mapped in the sublayer to layer mapping file
    valid_sublayers = sublayers_df['sublayer'].tolist()
    assert set(layers) - set(valid_sublayers) == set(), f'Invalid layer names, valid layer names are {valid_sublayers}'
    
    # Read in Table 1 and Table 2 data
    table1_filename = 'table1_data.csv'
    table2_filename = 'table2_data.csv'
    table1_fullpath = os.path.join(path_to_tables, table1_filename)
    table2_fullpath = os.path.join(path_to_tables, table2_filename)
    table1_df = pd.read_csv(table1_fullpath)
    table1_df = table1_df.set_index('stat') 
    table2_df = pd.read_csv(table2_fullpath)
   
    # Check that layer names match valid Table 2 layers
    assert set(sublayers_df['layer'].tolist()) - set(table2_df['layer'].tolist()) == set(), f'Invalid layer to Table 2 mapping'

    # Call generate_cell_groups to get df with neuron group names, numbers & positions by layer and type
    #Requested V1 area in mm2
    V1area = VFradius2V1area(requestedVFradius, center_ecc)

    # Get V1 proportion for simulation so that we can get N cells from the total N cells in V1 layers in Table 2
    V1total_area = table1_df.loc['mean','V1']
    v1_proportion = V1area / V1total_area

    # Get inhibitory cell type proportions from Markram rat somatosensory data
    cell_type_filename = 'layer_download.json'
    cell_type_fullpath = os.path.join(path_to_tables, 'hbp_data', cell_type_filename)
    cell_type_df = pd.read_json(cell_type_fullpath)
    inhibitory_proportions_markram_df = get_markram_cell_type_proportions(cell_type_df)

    # Get inhibitory cell type proportions from Allen institute data
    cell_type_allen_filename = 'sample_annotations.csv'
    cell_type_allen_fullpath = os.path.join(path_to_tables, 'allen_data', cell_type_allen_filename)
    cell_type_allen_df = pd.read_csv(cell_type_allen_fullpath)
    # Select appropriate classifying parameters from df
    cell_type_allen_relevant_columns_df = cell_type_allen_df[['region_label', 'class_label', 'subclass_label', 'cortical_layer_label']]
    # Get unique labels for comparison
    # valid_regions=set(cell_type_allen_relevant_columns_df['region_label']) # In case u want to compare regions later
    # valid_layers=set(cell_type_allen_relevant_columns_df['cortical_layer_label']) # In case u want to compare regions later
    # valid_class_labels=set(cell_type_allen_relevant_columns_df['class_label']) # In case u want to compare regions later
    # valid_subclass_labels=set(cell_type_allen_relevant_columns_df['subclass_label']) # In case u want to compare regions later
    # Select visual cortex according to region label
    cell_type_allen_V1_df = cell_type_allen_relevant_columns_df.loc[cell_type_allen_relevant_columns_df['region_label'] == 'V1C']
    cell_type_allen_V1clean_df = cell_type_allen_V1_df.drop(columns=['region_label'])
 
    # Excitatory proportions are not currently used
    inhibitory_proportions_allen_df, excitatory_proportions_allen_df = get_allen_cell_type_proportions(cell_type_allen_V1clean_df)

    # Get layer name mapping
    layer_name_mapping_filename = 'layer_name_mapping.xlsx'
    layer_name_mapping_fullpath = os.path.join(path_to_tables, layer_name_mapping_filename)
    layer_name_mapping_df = pd.read_excel(layer_name_mapping_fullpath)

    # Get df with neuron group names, numbers & positions by layer and type
    NG_df = generate_cell_groups(v1_proportion, layers, table2_df, inhibitory_proportions_allen_df, layer_name_mapping_df)

    

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

