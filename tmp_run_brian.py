from brian2 import *                                                                                                                                                                                

Nneurons=161
indices = arange(79,82,1) 
times = array([1, 2, 3])*ms                                                                                                                                                                         
G = SpikeGeneratorGroup(Nneurons, indices, times,period=5*ms)                                                                                                                                              

spike_mon = SpikeMonitor(G)                                                                                                                                                                         

run(1*second)                                                                                                                                                                                       
plt.scatter(spike_mon.t,spike_mon.i,s=1);plt.grid(which='major', axis='x');plt.show()                                                                                                               

