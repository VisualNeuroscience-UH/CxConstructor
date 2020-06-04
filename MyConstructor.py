from CxConstructor import *

'''
Start of user input
Copy and comment/uncomment examples/your own versions by need. If python gives exception, look first your own syntax below.
'''

area_name='V1' # Don't change this.
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
requested_layers=['L1', 'L23', 'L4A','L4B', 'L4CA','L5','L6'] # You should be able start from L4CA or B alone for testing
# requested_layers=['L1', 'L4C', 'L6'] # You should be able start from L4CA or B alone for testing

# Here are some examples
inhibitory_types = ['MC', 'BC']
# inhibitory_types = ['BC']
inhibitory_proportions={} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
# inhibitory_proportions = {  
# 'L1': [.5, .5], 
# 'L4C': [1, 0],
# 'L6': [.5, .5]}

# inhibitory_proportions = {  
# 'L1': [1, 0], 
# 'L23': [.5, .5], 
# 'L4A': [.5, .5], 
# 'L4B': [.5, .5], 
# 'L4C': [1, 0],
# 'L5': [.5, .5], 
# 'L6': [.5, .5]}

# Excitatory proportions are given by hand here. 
# The list length for each layer must match the N types and should sum to approximately 1.
# The layers must match the requested layers
excitatory_types = ['SS', 'PC1']
# excitatory_types = ['SS'] # IT is one of th Allen excitatory types. SS = spiny stellate and this should exist in all physiological files
# excitatory_proportions = {} # Leaving this empty will produce 1/N types proportions for each layer if cell_type_data_source = ''
# excitatory_proportions = {  
# 'L1': [1, 0],
# 'L4C': [1, 0],
# 'L6': [.1, .9]}
excitatory_proportions = {  
'L1': [1, 0], 
'L23': [0, 1], 
'L4A': [.5, .5], 
'L4B': [.5, .5], 
'L4CA': [1, 0],
'L5': [0, 1], 
'L6': [.1, .9]}

# Read in anat csv to start with. Check for config_files folder for valid inputs.
# If replace_existing_cell_groups flag is False, your groups will be added to current groups. This allows building system stepwise
replace_existing_cell_groups = True # 
anatomy_config_file_name = 'pytest_anatomy_config.csv'
physiology_config_file_name = 'pytest_physiology_config.csv' # anatomy and physiology filenames sometimes diverge

# Own additional data files
neuron_group_ephys_templates_filename = 'neuron_group_ephys_templates.xlsx'

# NOT FUNCTIONAL AT THE MOMENT 4.6.2020 #
# Activate correct selection. These data provide proportions of distinct cell groups in each layer.
# Markram data, rat S1
# cell_type_data_source = 'HBP'

# # Allen data, human V1
# cell_type_data_source = 'Allen'
# END OF NOT FUNCTIONAL AT THE MOMENT 4.6.2020 #

# No data source. Use this eg for single excitatory and single inhibitory types. 
cell_type_data_source = ''

request_monitors = '[Sp]'
n_background_inputs_for_excitatory_neurons = 630
n_background_inhibition_for_excitatory_neurons = 290    
n_background_inputs_for_inhibitory_neurons = 500
n_background_inhibition_for_inhibitory_neurons = 180    

# Connections. You can use only a subset of layers in ni csv data matching your requested layers or alternatively
# you can use all csv data, including the sublayers and CO compartments.
use_all_csv_data = True


###############################################
############## END OF USER INPUT ##############
###############################################

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

# Set Config class variables
Config.anatomy_config_df = read_config_file(os.path.join(PATH_TO_CONFIG_FILES,anatomy_config_file_name))
Config.physiology_config_df = read_config_file(os.path.join(PATH_TO_CONFIG_FILES,physiology_config_file_name))  
Config.replace_existing_cell_groups = replace_existing_cell_groups  
Config.table1_df = Config.get_neuroinformatics_data(TABLE1_DATA_FILENAME, set_index='stat')
Config.table2_df = Config.get_neuroinformatics_data(TABLE2_DATA_FILENAME, set_index='layer')
Config.use_all_csv_data = use_all_csv_data

V1 = Area(area_name=area_name, requestedVFradius=requestedVFradius, center_ecc=center_ecc, requested_layers=requested_layers)

# Add anatomy and physiology config files to start with
group_object = Groups(V1, requested_cell_types_and_proportions, cell_type_data_source, cell_type_data_folder_name, 
            cell_type_data_file_name, request_monitors,requested_background_input)

Conn_new = Connections(V1, group_object, use_all_csv_data)

# Write anatomy out
Config.write_config_files(Config.anatomy_config_df, anatomy_config_file_name, xlsx=True) # Original as excel file
Config.write_config_files(Conn_new.anatomy_config_df_new_groups_new_synapses, anatomy_config_file_name[:-4] + '_cxc', csv=True, xlsx=True)

# Write physiology out
Config.write_config_files(Config.physiology_config_df, physiology_config_file_name, xlsx=True) # Original as excel file
Config.write_config_files(group_object.physiology_df_with_subgroups, physiology_config_file_name[:-4] + '_cxc', csv=True, xlsx=True)
