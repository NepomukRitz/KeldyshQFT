import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex


#Load Vertex correction-Data:	
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Conductance/"
s_file="Precomputation_cond_vertex_cont.mat"

s=s_folder + s_file

freq_pre = np.real(ctn.load_as_np(s,"freq_pre"))[0,:]
gp = ctn.load_as_np(s,"gp");
gx = ctn.load_as_np(s,"gx");


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(freq_pre,np.real(gp[0,:,2,2]),color='blue',marker='o')
axis.plot(freq_pre,np.imag(gp[0,:,2,2]),color='red',marker='o')
axis.plot(freq_pre,np.real(gx[0,:,2,2]),color='green',marker='o')
axis.plot(freq_pre,np.imag(gx[0,:,2,2]),color='orange',marker='o')


plt.show()


