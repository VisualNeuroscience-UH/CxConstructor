import numpy as np
from matplotlib import pyplot as plt
from utilities import *

def newest(path):
    files = [f for f in os.listdir(path) if 'results' in f]
    paths = [os.path.join(path, basename) for basename in files]
    return max(paths, key=os.path.getctime)

latest_filename = newest(r'C:\Users\Simo\Laskenta\Models\SchwabeModel\build')
print(latest_filename)
data = getData(latest_filename)

t=data['vm_all']['NG1_SS_L2']['t']

vm1=data['vm_all']['NG1_SS_L2']['vm']

vm2=data['vm_all']['NG2_MC_L2']['vm']

vm3=data['vm_all']['NG3_SS_L6']['vm']

vm4=data['vm_all']['NG4_L2_SS_autoconn_L2']['vm']

plt.figure(1);plt.plot(t[2000:4000],vm1[2000:4000,0:-1:20]);

plt.figure(2);plt.plot(t[2000:4000],vm2[2000:4000,0:-1:20]);

plt.figure(3);plt.plot(t[2000:4000],vm3[2000:4000,0:-1:20]);

plt.figure(4);plt.plot(t[2000:4000],vm4[2000:4000,0:-1:20]);

plt.show()