import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

s="/gpfs/work/hmu26/hmu261/DATA/Unit_tests/Conductance_lead_system/cond_lead_sys.mat"

#vertex_contribution = ctn.load_as_np(s,"vertex_contribution");
freq_subst = ctn.load_as_np(s,"freq_subst");
freq = ctn.load_as_np(s,"freq");
Integrated_vertex_contribution = ctn.load_as_np(s,"Integrated_vertex_contribution");
Old_integrated_vertex_contribution = ctn.load_as_np(s,"Old_integrated_vertex_contribution");

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))

#axis[0].plot(freq_subst[0,:],np.real(vertex_contribution[0,:]), color='blue', marker='o')
#axis[0].plot(freq[0,:],np.real(vertex_contribution[0,:]), color='blue', marker='o')
axis[0].plot(freq[0,:],np.real(Old_integrated_vertex_contribution[0,:,15,15]), color='blue', marker='o')
axis[0].plot(freq[0,:],np.real(Integrated_vertex_contribution[0,:,15,15]), color='red', marker='o')

plt.show()

