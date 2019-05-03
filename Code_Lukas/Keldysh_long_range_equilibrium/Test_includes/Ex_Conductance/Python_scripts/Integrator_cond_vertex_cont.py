import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex


#Load Vertex correction-Data:	
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Conductance/"
s_file="Integrator_cond_vertex_cont.mat"

s=s_folder + s_file

Omega = np.real(ctn.load_as_np(s,"Omega"))[0,:]
N_omega = np.real(ctn.load_as_np(s,"N_omega"))[0,0].astype(int)
wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
N = np.real(ctn.load_as_np(s,"N"))[0,0].astype(int)
Vertex_contribution = ctn.load_as_np(s,"Vertex_contribution")
A = ctn.load_as_np(s,"A")

print(A[0,0,:,:])

print("N=",N)
print("Vertex_contribution=", Vertex_contribution[0,1,N,N]); 


#fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#axis.plot(Omega,np.real(Vertex_contribution[0,:,N,N]),color='blue',marker='o')
#axis.plot(Omega,np.imag(Vertex_contribution[0,:,N,N]),color='red',marker='o')
#
#
#plt.show()


