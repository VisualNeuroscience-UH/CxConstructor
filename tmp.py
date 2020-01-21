import numpy as np
from matplotlib import pyplot as plt
from utilities import *

# saveSchwabePositions('model06positions',dir_name=r'C:\Users\Simo\Laskenta\Models\SchwabeModel\positions')

# data=getData(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\positions\model06positions.gz')
data=getData(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build\build_connections_20200106_1330045.gz')

# Load connections with correct lambda
# data = getData(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build\lambda_connections.gz')

# # Load connections with correct anatomical model
# data = getData(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build\model_connections.gz')

# # Load results 
# data = getData(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build\build_results_20200106_0609316_nan500.gz')

data = buildSchwabeConnections(data)

saveSchwabeData(data, 'model_connections_200106',dir_name=r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build')

# # ng1= data['positions_all']['w_coord']['NG1_SS_L4']
# # ng0= data['positions_all']['w_coord']['NG0_relay_vpm']
# # ng2= data['positions_all']['w_coord']['NG2_BC_L4']
# # ng3= data['positions_all']['w_coord']['NG3_SS_L2']

# # ng0a=np.real(np.asarray(ng0))
# # ng1a=np.real(np.asarray(ng1))
# # ng2a=np.real(np.asarray(ng2))
# # ng3a=np.real(np.asarray(ng3))

# # plt.figure()
# # plt.scatter([ng0a],[np.ones(ng0a.shape)*0])
# # plt.scatter([ng1a],[np.ones(ng1a.shape)*1])
# # plt.scatter([ng2a],[np.ones(ng2a.shape)*2])
# # plt.scatter([ng3a],[np.ones(ng3a.shape)*3])
# # plt.show()
