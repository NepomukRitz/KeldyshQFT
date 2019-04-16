import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys
import One_particle_conductance_to_numpy as Opctn 

s = "/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Conductance/Compare_vertex_correction.mat"

vertex_correction_list = ctn.load_as_np(s,"vertex_correction_list")
VertCorr_list = ctn.load_as_np(s,"VertCorr_list")
external_freq_list = ctn.load_as_np(s,"external_freq_list")[0,:]

print("np.amax(np.absolute(vertex_correction_list - VertCorr_list))=",np.amax(np.absolute(vertex_correction_list - 2.0*VertCorr_list)))
fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#axis.plot(mu_val,cond_val_up+cond_val_down,color='blue',marker='o')
#axis.plot(mu_val_NL10nle,cond_val_up_NL10nle+cond_val_down_NL10nle,color='red',marker='o')
#axis.plot(mu_val_NL2,cond_val_up_NL2+cond_val_down_NL2,color='green',marker='o')
#axis.plot(mu_val_NL5,cond_val_up_NL5+cond_val_down_NL5,color='red',marker='o')
#axis.plot(mu_val_NL10,cond_val_up_NL10+cond_val_down_NL10,color='orange',marker='o')

axis.plot(external_freq_list,np.real(vertex_correction_list[0,:,2,1]),color='blue',marker='o')
axis.plot(external_freq_list,np.imag(vertex_correction_list[0,:,2,1]),color='red',marker='o')
axis.plot(external_freq_list, np.real(VertCorr_list[0,:,2,1]),color='green',marker='o')
axis.plot(external_freq_list, np.imag(VertCorr_list[0,:,2,1]),color='orange',marker='o')

#axis.set_xlim([-5,5])

plt.show()



