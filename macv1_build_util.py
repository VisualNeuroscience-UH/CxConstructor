import numpy as np
from matplotlib import pyplot as plt
import os
import pandas as pd
import pdb

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
def generate_cell_groups(requestedV1area, center_ecc, layers):
    pass
'''
generate_cell_groups function
-for each layer in macaque V1, set the neuron groups, types and numbers
    -unit N / mm2 / layer
    -requires Table2, cellTypeProportions, V1areaTotal
    -input requestedV1area, center_ecc, layers
    -set the V1 area in mm2 (total) according to M function
    -read in the requested V1 area
    -read in Table 2 as df (N neurons/layer, N synapses / layer, %inhibitory/layer)
    -read in estimates of cell type proportions (by layer) from separate file
    -map the Table2 data into cell groups
    -calculate density for each group
    -set positions by layer and type
    -return df with neuron group names, numbers & positions by layer and type 
'''

def cell_group_row():
    pass
'''
cell_group_row function
We want to keep this separate method, because the syntax of the output csv file may change
-return row_type,idx,number_of_neurons,neuron_type,neuron_subtype,layer_idx,net_center,monitors,n_background_inputs,n_background_inhibition
'''

def generate_synapses():
    pass
'''
generate_synapses function
-for each origin and target group, set receptor,pre_syn_idx,post_syn_idx by layer and compartment,syn_type,p,n,
-return df with receptor,pre_syn_idx,post_syn_idx,syn_type,p,n, 
'''

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

def ni2anat(layers, relative_depth):
    pass
'''
ni2anat function
A filter for turning neuroinfo csv data into connections in anatomical csv file
-read ni csv into dataframe
-input layers,  
    
The sublayer to layer mapping file contains eg sublayers 4Cm, L6A etc. Cell proportions in sublayers are simply
1/(N sub layers). The blob/interblob relation is estimated to be 0.4/0.6, according to surface areas at 5 deg ecc (Livingstone 1984 JNeurosci)

-check their names for validity from the sublayer to layer mapping file
-input relative_depth (proportion cells) for sublayers
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

