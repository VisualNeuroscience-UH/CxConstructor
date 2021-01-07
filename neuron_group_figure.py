import os
import matplotlib.pyplot as plt 
import numpy as np
import utilities as ut

TÄHÄN JÄIT. TEE KUVA JOSSA ROWS ON NG JA COLS ON PARAMS

class showFigure:

    figure_size = (4,6)
    path_to_data = '/opt3/CxSimoWorkspace/Grossberg'
    filename = 'Grossberg_results_20200630_1710246.gz'

    def __init__(self, vm=True, raster=True, spectrum=True, coherence=True, transfer_entropy=True):

        # Set visible columns
        self.vm = vm
        self.raster = raster
        self.spectrum = spectrum
        self.coherence = coherence
        self.transfer_entropy = transfer_entropy
        self.ncols = vm + raster + spectrum + coherence + transfer_entropy

        # get data
        self.data = ut.getData(os.path.join(showFigure.path_to_data, showFigure.filename))

        # get visible rows
        self.neuron_groups = list(self.data['spikes_all'].keys())
        self.anat_config = self.data['Anatomy_configuration']
        self.nrows = len(self.neuron_groups)
        

    def passer(self):
        pass

    def myData(self, pathname, filename):
        data = ut.getData(os.path.join(pathname, filename))
        return data

def main():
    myFig = showFigure()
    
    # print(type(myFig.anat_config))
    print(myFig.ncols, myFig.nrows)

if __name__ == "__main__":
    main()