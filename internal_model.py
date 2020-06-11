'''
Build images of receptive fields and internal representations of simulated data.
'''

import numpy as np
import os
from utilities import getData
import pdb

def get_rf():

    # From connections, get neuron positions and weights from input group to target group
    # In case of second order neurons, get the first-level RF first, then use them to get the
    # second level weights. 

def get_internal_image(ng,simulation_filename):
    # Pick neuron group, get rf for each neuron, weigh the rf with response rate from simulation. 
    # Show internal image

    data = getData(simulation_filename)

    pdb.set_trace()