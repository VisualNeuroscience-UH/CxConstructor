# -*- coding: utf-8 -*-
'''
Construct macaque visual cortex model.

This module input your definitions after the if __name__ == "__main__":
and constructs model of macaque V1. It makes a query to neuroinfo csv data 
and tables from Vanni et al Cerebral Cortex 2020. It uses these data to create
cell groups and connections into anatomical csv file for CxSystem2 package.

Classes
    Config: This class instance contains general objects, concerning all areas. It is here for taking data around and including general methods.
    Area:   This class contains area-level data and methods
    Group:  Neuron group level data and methods
    Connections:    Generate synapses object, which includes the new anatomy df with connections

Design
    -read ni csv into dataframe, Table1 and Table2
    -input layers,  
        
    The layer_name_mapping excel file contains eg sublayers 4Cm, L6A etc. Cell proportions in sublayers are simply
    1/(N sub layers). The blob/interblob relation is estimated to be 0.4/0.6, according to surface areas at 5 deg ecc (Livingstone 1984 JNeurosci)

    -check requested names for validity

    -read in Table 1 and Table 2 data
    -map the N layers to valid Table 2 layers, 

    -read in anat csv to start with

    -call generate_cell_groups to get df with neuron group names, numbers by layer and type/subtype
    -call cell_group_row for anat csv row output for cell groups

    -map each row in the ni csv file into layer of origin and termination, set strength
    -set origin groups according to subgroups
    -set target groups according to subgroups
    -call generate_synapses function
    -call synapse_row

    -write cell groups and synapses into the existing anat csv file with _cxc suffix

    -call show_cortex

'''
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

# TODO Instead of passing objects from class to class, consider inheritance



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
        layer_name_mapping_df_orig = CxC.read_data_from_tables(path_to_tables, 'layer_name_mapping_V1.xlsx')
        # Drop comments
        layer_name_mapping_df_orig = layer_name_mapping_df_orig.drop(columns='comments')
        self.layer_name_mapping_df_orig = layer_name_mapping_df_orig
        self.PC_apical_dendrites = CxC.read_data_from_tables(path_to_tables, 'PC_apical_dendrites.xlsx')

        # Check layer names for validity: are they mapped in the sublayer to layer mapping file
        # valid_layers = layer_name_mapping_df_orig['allowed_requested_layers'].tolist()
        layer_name_mapping_df_orig_search_columns = ['down_mapping1', 'down_mapping2', 'down_mapping3', 'csv_layers']
        self.layer_name_mapping_df_orig_search_columns = layer_name_mapping_df_orig_search_columns

        valid_layers = self.get_valid_layers(   layer_name_mapping_df_orig, 
                                                layer_name_mapping_df_orig_search_columns)
        assert set(requested_layers) - set(valid_layers) == set(), f'Invalid layer names, valid layer names are {valid_layers}'
        self.valid_layers = valid_layers

        # Map requested_layers to layer_name_mapping_df_groups. Calculate proportions for groups  
        self.layer_name_mapping_df_groups, self.layer_name_mapping_df_full = \
                                        self.map_requested_layers2valid_layers(
                                        requested_layers, layer_name_mapping_df_orig, 
                                        layer_name_mapping_df_orig_search_columns)

        if area_name=='V1':
            # Get proportion V1 of total V1
            #Requested V1 area in mm2
            V1area = self._VFradius2V1area(requestedVFradius, center_ecc)

            # Get V1 proportion for simulation so that we can get N cells from the total N cells in V1 layers in Table 2
            V1total_area = CxC.table1_df.loc['mean','V1']
            self.area_proportion = V1area / V1total_area
        else:
            self.area_proportion = 1 # Not implemented yet for other areas

    def map_requested_layers2valid_layers(self, requested_layers, layer_name_mapping_df_orig, search_columns):
        '''
        Map requested_layers to layer_name_mapping_df_groups. Calculate proportions for groups.
        Return df with one layer per requested layer with proportions calculated
        '''
        # Calculate proportions
        # Loop through requested layers
        # Loop through layer_name_mapping_df_orig
        # Is requested layer in down_mapping1, down_mapping2, down_mapping3 or in csv_layers columns?
        # Check if this csv_layer is already accounted for: 
        #   is any down mapping column names already in picked csv list?
        # If yes, discard this layer and continue
        # If no, append corresponding rows to layer_name_mapping_df
        # mark this csv_layer picked for layer_name_mapping_df
        layer_name_mapping_df = pd.DataFrame(columns=layer_name_mapping_df_orig.columns)

        # Add columns for layer idx and requested layers
        layer_name_mapping_df = layer_name_mapping_df.join(pd.DataFrame(columns=['layer_idx']))
        layer_name_mapping_df = layer_name_mapping_df.join(pd.DataFrame(columns=['requested_layers']))
        layer_name_mapping_df_full = layer_name_mapping_df

        # Loop through requested layers
        for current_layer_requested_idx, current_layer_requested in enumerate(requested_layers):

            # Loop through layer_name_mapping_df_orig
            for index, row in layer_name_mapping_df_orig.iterrows():
                
                # Is requested layer in this row's search_columns?
                if current_layer_requested in row[search_columns].values:
                    
                    # Check if current csv_layer is already accounted for: 
                    #   are any of the down mapping column names (row) of the current row 
                    #   already in the csv_layers column of the layer_name_mapping_df?

                    #   if set does not change (no matching names in the csv_layers column), 
                    #   append row to layer_name_mapping_df, otherwise continue
                    if (set(row[search_columns].values) - set(layer_name_mapping_df['csv_layers'].values)) == set(row[search_columns].values):
                        layer_name_mapping_df = layer_name_mapping_df.append(row)
                        # Add layer_idx to last row
                        layer_name_mapping_df.loc[index,'layer_idx'] = current_layer_requested_idx + 1
                        layer_name_mapping_df.loc[index,'requested_layers'] = current_layer_requested
                
                    # All layers are included in full for full csv mapping
                    layer_name_mapping_df_full = layer_name_mapping_df_full.append(row)
                    layer_name_mapping_df_full.loc[index,'layer_idx'] = current_layer_requested_idx + 1
                    layer_name_mapping_df_full.loc[index,'requested_layers'] = current_layer_requested

        # Set layer_idx to integer type
        layer_name_mapping_df['layer_idx'] = layer_name_mapping_df['layer_idx'].astype('int32')
        layer_name_mapping_df_full['layer_idx'] = layer_name_mapping_df_full['layer_idx'].astype('int32')

        # # Add proportions to layer_name_mapping_df_full
        # # layer_name_mapping_df_full = layer_name_mapping_df_full.join(pd.DataFrame(columns=['sub_proportion']))

        # layer_name_mapping_df_full['sub_proportion'] = layer_name_mapping_df_orig['sub_proportion']

        # Check that sums of cells do not exceed table2 counts
        # assert layer_name_mapping_df['sub_proportion'].sum(axis=0) <= current_layer_requested_idx + 1, 'Sums of cells will exceed table2 counts'

        return layer_name_mapping_df, layer_name_mapping_df_full

    def get_valid_layers(self, layer_name_mapping_df_orig, columns_for_search):
        '''
        Search layer_name_mapping_V1 through columns csv_layers, and down_mapping1-3.
        Any unique name is valid for request
        '''
        # Columns for search
        # columns_for_search = ['down_mapping1', 'down_mapping2', 'down_mapping3', 'csv_layers']
        valid_layers = pd.unique(layer_name_mapping_df_orig[columns_for_search].values.ravel('K'))
        return valid_layers

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
    '''
    Neuron group level data and methods
    '''
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

        # Get proportions of inhibitory and excitatory neurons in each layer
        # Valid EIflag 'Glutamatergic' and 'GABAergic'
        self.inhibitory_proportions_df = self.get_proportions_df(   'GABAergic',inhibitory_proportions, inhibitory_types, requested_layers, 
                                                                    cell_type_data_source, cell_type_data_folder_name, cell_type_data_file_name)
        
        self.excitatory_proportions_df = self.get_proportions_df(   'Glutamatergic',excitatory_proportions, excitatory_types, requested_layers, 
                                                                    cell_type_data_source, cell_type_data_folder_name, cell_type_data_file_name)
        
        # Map cell groups to requested layers
        # Choose layer mappings according to requested layers. 
        # TODO If an entry has two values separated by ";", the two values must be averaged

        # # layer_name_mapping_df_orig = area_object.layer_name_mapping_df_orig # Should be acquired from Area class instance object
        # # # Map requested layers to allowed neuroinformatics sublayers.
        # # layer_mapping_df = layer_name_mapping_df_orig.loc[layer_name_mapping_df_orig['allowed_requested_layers'].isin(requested_layers)]
        # # # Add numerical layer index to layer_mapping_df
        # # layer_idxs = np.arange(len(requested_layers)) + 1
        # # layer_mapping_df['layer_idx'] = layer_idxs.tolist()
        self.layer_mapping_df = area_object.layer_name_mapping_df_groups

        # # self.layer_mapping_df = layer_mapping_df

        # Get df with neuron groups for anatomy df, return new df to object
        self.anatomy_config_df_new = self.generate_cell_groups(CxC, area_object, requested_cell_types_and_proportions)

    def generate_cell_groups(self, CxC, area_object, requested_cell_types_and_proportions):
        '''
        Generate_cell_groups function
        For each layer in macaque V1, set the neuron groups, types and numbers
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

        # Count number of new groups, ie number of new rows
        index_excitatory_s = excitatory_proportions_df.fillna(0).astype(bool).sum()
        index_inhibitory_s = inhibitory_proportions_df.fillna(0).astype(bool).sum()
        N_rows =  index_excitatory_s.sum() + index_inhibitory_s.sum()
        
        # Generate df for holding the neuron groups. Add indices here to preallocate memory
        indices = start_index + np.arange(N_rows)      
        NG_df = pd.DataFrame(columns=cell_group_columns, index=indices)
        
        # Now we go for the cell groups. First let's set the identical values
        NG_df.row_type = 'G'
        NG_df.idx = np.arange(start_cell_group_index, start_cell_group_index + N_rows)
        NG_df.net_center = NG_df.noise_sigma = NG_df.gemean = NG_df.gestd = NG_df.gimean = NG_df.gistd = '--'
        NG_df.monitors = monitors

        # Add one layer at a time starting from L1
        current_group_index = start_index

        for layer in requested_layers:
            # current_layer_proportion_inhibitory = table2_df.loc[layer]['percent_inhibitory'] / 100
            current_layer_proportion_inhibitory = self.calc_proportion_inhibitory(layer, table2_df, layer_mapping_df)
            # current_layer_N_neurons_area = table2_df.loc[layer]['n_neurons_10e6'] * 1e6
            current_layer_N_neurons_area = self.calc_N_neurons(layer, table2_df, layer_mapping_df)
            current_layer_N_excitatory_neurons = (1 - current_layer_proportion_inhibitory) * area_proportion * current_layer_N_neurons_area
            current_layer_N_inhibitory_neurons = current_layer_proportion_inhibitory * area_proportion * current_layer_N_neurons_area
            # layer_idx = layer_mapping_df.loc[layer_mapping_df['allowed_requested_layers']==layer]['layer_idx'].values
            layer_idx = layer_mapping_df.loc[layer_mapping_df['requested_layers'] == layer, 'layer_idx'].values[0]

            # Add excitatory cell groups
            for current_group in excitatory_proportions_df.index.values:
                # Jump to next group if proportion is 0
                if excitatory_proportions_df.loc[current_group][layer]==0:
                    continue
                # go through all columns which require individual values. 
                # TODO Add from connection_csv_sublayer_to_table2_mapping.csv subproportion here to number of neurons
                NG_df.loc[current_group_index,'number_of_neurons'] =  np.round(excitatory_proportions_df.loc[current_group][layer] * current_layer_N_excitatory_neurons)
                NG_df.loc[current_group_index,'neuron_type'] = current_group
                # Use layers as neuron subtype prefixes
                NG_df.loc[current_group_index,'neuron_subtype'] = layer + '_' + current_group
                # For PC, map apical dendrite extent from table to layer index
                if current_group == 'PC':
                    # TODO: NOT IMPLEMENTED YET
                    pass
                    # NG_df.loc[current_group_index,layer_idx] = layer_mapping_df
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

    def calc_N_neurons(self, current_layer, table2_df, layer_mapping_df):
        # When more or less than than one layer maps to requested layer, one cannot directly query table 2
        # Instead you need to sum according to sub_proportion
        # from layer_mapping_df, choose current_layer rows in requested_layer_column
        # TODO Less than one layer not implemented yet

        table2to_sub_proportions_df = layer_mapping_df.loc[layer_mapping_df['requested_layers'] == 
                                        current_layer, ['table2_df', 'sub_proportion']]

        # Calculate sum, 
        # sigma_i(sub_proportion_i * n_neurons_i) * 1e6

        sum_neurons = 0
        for idx, table2_requested_layer in enumerate(table2to_sub_proportions_df['table2_df'].values):
            sum_neurons += table2to_sub_proportions_df.loc[:,'sub_proportion'].values[idx] * \
                                table2_df.loc[table2_requested_layer,'n_neurons_10e6'] 

        # Scale
        sum_neurons_scaled = sum_neurons * 1e6

        return sum_neurons_scaled

    def calc_proportion_inhibitory(self, current_layer, table2_df, layer_mapping_df):
        # When more or less than than one layer maps to requested layer, one cannot directly query table 2
        # Instead you need to average according to sub_proportion
        # from layer_mapping_df, choose current_layer rows in requested_layer_column
        # TODO Less than one layer not implemented yet
        # TODO remove quickfix where 4CM is accounted by 4CB in down_map3 of table layer_name_mapping_V1

        table2to_sub_proportions_df = layer_mapping_df.loc[layer_mapping_df['requested_layers'] == 
                                        current_layer, ['table2_df', 'sub_proportion']]

        # Calculate weighed average, 
        # sigma_i(sub_proportion_i * n_neurons_i * percent_inhibitory_i)/sigma_i(sub_proportion_i * n_neurons_i)

        useful_numerator = 0
        useful_denominator = 0
        for idx, table2_requested_layer in enumerate(table2to_sub_proportions_df['table2_df'].values):
            useful_numerator += table2to_sub_proportions_df.loc[:,'sub_proportion'].values[idx] * \
                                table2_df.loc[table2_requested_layer,'n_neurons_10e6'] * \
                                table2_df.loc[table2_requested_layer,'percent_inhibitory']
            useful_denominator +=   table2to_sub_proportions_df.loc[:,'sub_proportion'].values[idx] * \
                                    table2_df.loc[table2_requested_layer,'n_neurons_10e6']
        percentage_inhibitory = useful_numerator / useful_denominator

        # From percentage to proportion
        proportion_inhibitory = percentage_inhibitory / 100

        return proportion_inhibitory

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
    '''
    Generate synapses object, which includes the new anatomy df with connections
    '''

    def __init__(self, CxC, area_object, NG_new, use_all_csv_data):
        
        # Set parameters to connections object
        self.replace_existing_cell_groups = CxC.replace_existing_cell_groups
        
        # Read data from files.
        # Read ni csv into dataframe
        exc_df = CxC.read_data_from_tables(path_to_ni_csv, 'connections_local_excitatory.csv')
        inh_df = CxC.read_data_from_tables(path_to_ni_csv, 'connections_local_inhibitory.csv')
        area_name = area_object.area_name
        requested_layers = area_object.requested_layers
        # layer_mapping_df = NG_new.layer_mapping_df
        # layer_name_mapping_df_orig = area_object.layer_name_mapping_df_orig
        layer_name_mapping_df_groups = area_object.layer_name_mapping_df_groups
        layer_name_mapping_df_full = area_object.layer_name_mapping_df_full


        self.anatomy_config_df_new_groups = NG_new.anatomy_config_df_new
        
        # Map exc_df and inh_df to valid format inhibitory and excitatory connections in each layer
        self.excitatory_connections_df = self.get_local_connection_df(exc_df, area_name, 
                                            layer_name_mapping_df_full, requested_layers,
                                            layer_name_mapping_df_groups, use_all_csv_data)

        self.inhibitory_connections_df = self.get_local_connection_df(inh_df, area_name, 
                                            layer_name_mapping_df_full, requested_layers,
                                            layer_name_mapping_df_groups, use_all_csv_data)
        pdb.set_trace()
        # Generate connections
        self.anatomy_config_df_new_groups_new_synapses = self.generate_synapses()


    def generate_synapses(self):
        '''
        generate_synapses function
        -for each origin and target group, set receptor,pre_syn_idx,post_syn_idx by layer and compartment,syn_type,p,n,
        -return df with receptor,pre_syn_idx,post_syn_idx,syn_type,p,n, 

        '''

        # Unpack for this method
        replace_existing_cell_groups = self.replace_existing_cell_groups
        excitatory_connections_df = self.excitatory_connections_df 
        inhibitory_connections_df = self.inhibitory_connections_df 
        anatomy_config_df_new_groups = self.anatomy_config_df_new_groups
        # layer_mapping_df = self.layer_mapping_df

        # Get and set column names for connections
        connection_columns = anatomy_config_df_new_groups.loc[anatomy_config_df_new_groups.groupby([0]).get_group('S').index[0]-1,:].values
        existing_connection = anatomy_config_df_new_groups.groupby([0]).get_group('S')
        existing_connection.columns = connection_columns

        # Get starting index for connections and row indices for anat df
        if CxC.replace_existing_cell_groups:
            # start_connection_index = 1 # Reserve 0 for input group
            start_index = anatomy_config_df_new_groups.groupby([0]).get_group('S').index[0]
        else:
            # start_connection_index = int(existing_connection['idx'].values[-1]) + 1
            start_index = anatomy_config_df_new_groups.groupby([0]).get_group('S').index[-1] + 1

        # Count number of new connections, ie number of new rows

        pdb.set_trace()

        pass

    def get_local_connection_df(self, ni_df, area_name, layer_name_mapping_df_full, 
                                requested_layers,layer_name_mapping_df_groups,
                                use_all_csv_data):
        '''
        Turn neuroinformatics connection df to useful connection_df

        Select area
        Turn layer names in csv into requested layer names

        Three categories of csv connections. 
        1) requested layer names found: Use directly csv layer name
        2) requested layer names not found: Down-map from requested layer to csv         
        3) do not use if flag use_all_csv_data = False

        Create map from D, M, S to connection probabilities. Allows use of proportions
        Separate proportions for the pre- and postsynaptic side, multiply for final connections

        '''

        # Select area
        connection_df_ni_names = ni_df.groupby(['FromArea']).get_group(area_name)

        # Drop unnecessary area columns
        connection_df_ni_names = connection_df_ni_names.drop(columns=['FromArea', 'ToArea'])

        # Strip references
        connection_df_ni_names = connection_df_ni_names.drop(columns=['References'])
        
        # Duplicate csv layer names for proportion columns
        connection_df_ni_names['FL_proportions'] = connection_df_ni_names['FromLayer']
        connection_df_ni_names['TL_proportions'] = connection_df_ni_names['ToLayer']

        if use_all_csv_data:
            # layer_name_mapping_df_groups maps requested layer names, layer idx and csv layer names
            csv2idx_dict = dict(zip(layer_name_mapping_df_full.csv_layers, 
                                    layer_name_mapping_df_full.layer_idx))
            csv2proportion_dict = dict(zip(layer_name_mapping_df_full.csv_layers, 
                                    layer_name_mapping_df_full.sub_proportion))
        else:
            csv2idx_dict = dict(zip(layer_name_mapping_df_groups.csv_layers, 
                                    layer_name_mapping_df_groups.layer_idx))
            csv2proportion_dict = dict(zip(layer_name_mapping_df_groups.csv_layers, 
                                    layer_name_mapping_df_groups.sub_proportion))

        # Replace csv layer names with idx. The FL_proportions and TL_proportions indicate
        # relative layer thickness in comparison to Table 2 layers.
        connection_df_ni_names['FromLayer'] = connection_df_ni_names['FromLayer'].replace(csv2idx_dict)
        connection_df_ni_names['ToLayer'] = connection_df_ni_names['ToLayer'].replace(csv2idx_dict)
        connection_df_ni_names['FL_proportions'] = \
            connection_df_ni_names['FL_proportions'].replace(csv2proportion_dict)
        connection_df_ni_names['TL_proportions'] = \
            connection_df_ni_names['TL_proportions'].replace(csv2proportion_dict)

        # Now we just skip the rows not matching either pre- or postsynaptic layers
        layer_idxs = layer_name_mapping_df_groups.layer_idx.unique()
        connection_df_ni_names = connection_df_ni_names.loc[connection_df_ni_names['FromLayer'].isin(layer_idxs)]
        connection_df_ni_names = connection_df_ni_names.loc[connection_df_ni_names['ToLayer'].isin(layer_idxs)]


        '''
        ***Mapping from connection strength to connection probability***

        Create map from D, M, S to connection probabilities. Allows use of proportions
        In brian2 the synaptic p (probability of connection) parameter provides the 
        N connections / (NcellsG1 * NcellsG2), ie the number of connections from all possible pairs.
        The brian2 n-parameter multiplies each existing connection with integer number.
        
        The connection strength is the proportion of axonal sprouting/labelled cells in particular layer.
        We map connection strength value to numerical weight:
        D, dominant = 0.5 - 1, mean 0.75
        M, median = 0.1 - 0.5, mean 0.3
        S, sparse = eps - 0.1, mean 0.05        
        
        Local connections
        Intracellular tagging, shows single neuron structure. 
        We get:
        p = sw * (lw/0.5) * tw * apc
        p: probability of connections between any two cells (brian parameter). 
        sw: source group weight. N efferent synapses for this neuron group / mean N efferent synapses in cortical neurons
        lw: layer weight. The connection strength parameter from ni paper. The proportion of axonal sprouting in particular layer
        tw: target group weight. N presynaptic terminals from this neuron group / mean N presynaptic terminals from any neuron group
        tw 0 means avoidance of the postsynaptic group. tw 1 means Peter's rule
        apc: average probability of connection between two neurons, when efferent axon and afferent cell soma are in the same layer

        N pre and postsynaptic cells is already accounted for in NeuronGroups.
        sw, lw/0.5 and tw are all mean 1

        Because we do not know sw or tw for our neuron groups, this becomes
        p = (lw/0.5) * apc


       ***Weighing sublayers from ni cvs to probability***

        When we have multiple sublayers, we need to weigh the p value according estimated relative 
        sublayer thickness which we assume relates to the proportion of neurons in the sublayer from 
        all neurons in the layer. Eg for source sublayer1 to target sublayer1 connection , the p_partial value is
        p1 =  p_s1_t1 * sst1 * tst1
        p_s1_t1: probability of connection between source and target sublayers
        sst1: source sublayer thickness
        tst1: target sublayer thickness
        
        
        ***Scaling total p value according to included total sublayer thicknesses***
        
        The total p value for connecting two neuron groups becomes
        p_total = p1 + p2 + ... pN, where the partial p values indicate connection probabilities between each sublayer.
 
        Finally, include scaling by total proportions over all included sublayers. Full layer gives one and
        if u add sublayers A and B they include 0.5 each and BL and IB again 0.4 and 0.6 which sums to over 1.

        source scaling factor ssf = sst1 + sst2 + ... + sstN
        target scaling factor tsf = tst1 + tst2 + ... + tstN
        
        Final p value:

        p = p_total * 2 / (ssf + tsf)
        '''

        # Calculate partial p values for each sublayer to sublayer connection
        cw2p_dict = dict(zip(['D', 'M', 'S'],[0.75, 0.3, 0.05]))
        apc = 0.1
        sw = tw = 1
        connection_df_ni_names['Strength'] = connection_df_ni_names['Strength'].replace(cw2p_dict)
        connection_df_ni_names['p_partial'] = \
            sw * (connection_df_ni_names['Strength'] / 0.5) * tw * apc * \
            connection_df_ni_names['FL_proportions'] * connection_df_ni_names['TL_proportions']
            
        # Scale each unique set of partial connections to total FL_portion + TL_portion = 2
        # Change appropriate columns' data types to float. Mixed integers seemed to make these "objects"
        connection_df_ni_names = connection_df_ni_names.astype({'FL_proportions': 'float', 'TL_proportions': 'float', 'p_partial': 'float'})        
        
        # Group according to unique connections and sum all float columns
        connection_df_ni_names_unique = connection_df_ni_names.groupby(['FromLayer','ToLayer']).sum()
        connection_df_ni_names_unique = connection_df_ni_names_unique.reset_index()

        # I am not a perfectionist but I work with a mathematician and want to be exact, thus
        connection_df_ni_names_unique = connection_df_ni_names_unique.rename(columns={'p_partial':'p_total'})
        
        # Scale to final p value
        connection_df_ni_names_unique['p'] = (connection_df_ni_names_unique['p_total'] * 2) / \
            (connection_df_ni_names_unique['FL_proportions'] + connection_df_ni_names_unique['TL_proportions'])
        

        return connection_df_ni_names_unique


def show_cortex():
    '''
    show_cortex
    Visualize the neurons (and connections?) in 3D
    visimpl?
    '''    
    pass


if __name__ == "__main__":
    '''
    Start of user input
    Copy and comment/uncomment examples/your own versions by need. If python gives exception, look first your own syntax below.
    '''

    area_name='V1' # Don't change this.
    # requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA', 'L4CB','L5','L6'] # You should be able start from L4CA or B alone for testing
    requested_layers=['L1', 'L23', 'L4A','L4B', 'L4C','L5','L6'] # You should be able start from L4CA or B alone for testing
    # requested_layers=['L23', 'L4CA', 'L5','L6'] # You should be able start from L4CA or B alone for testing
    # requested_layers=['SG', 'L4', 'IG'] # You should be able start from L4CA or B alone for testing
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
    inhibitory_types = ['MC', 'BC']
    # inhibitory_types = ['BC']
    # inhibitory_proportions={} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
    inhibitory_proportions = {  
    'L1': [1, 0], 
    'L23': [.5, .5], 
    'L4A': [.5, .5], 
    'L4B': [.5, .5], 
    'L4C': [1, 0],
    'L5': [.5, .5], 
    'L6': [.5, .5]}
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
    # excitatory_types = ['SS', 'PC1', 'PC2']
    excitatory_types = ['SS'] # IT is one of th Allen excitatory types. SS = spiny stellate and this should exist in all physiological files
    excitatory_proportions = {} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
    # excitatory_proportions = {  
    # 'L1': [1, 0, 0], 
    # 'L23': [0, .5, .5], 
    # 'L4A': [.33, .33, .33], 
    # 'L4B': [.33, .33, .33], 
    # 'L4CA': [1, 0, 0],
    # 'L4CB': [1, 0, 0],
    # 'L5': [0, .5, .5], 
    # 'L6': [0, .5, .5]}

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

    # Connections
    use_all_csv_data = False

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
    # For debugging we write this to excel, too
    NG_new.anatomy_config_df_new.to_excel(os.path.join(path_to_config_files,anatomy_config_file_name[:-4] + suffix_for_new_file + '.xlsx'), header=False, index=False)
    # For cxsystem we write this to json, too
    NG_new.anatomy_config_df_new.to_json(os.path.join(path_to_config_files,anatomy_config_file_name[:-4] + suffix_for_new_file + '.json'))

    Conn_new = Connections(CxC, V1, NG_new, use_all_csv_data)

