import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys

s="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp/ex_dsfRG_L_0_Lu_0_N_30_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_0_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.000000_U0_0.650000_U1_0.000000_Xi_5.000000.mat"

U = np.real(ctn.load_as_np(s,"U"))
print(np.diag(U))
Lambda_initial = np.real(ctn.load_as_np(s,"Lambda_initial"))
Lambda_final = np.real(ctn.load_as_np(s,"Lambda_final"))
print(Lambda_initial)
print(Lambda_final)


#fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#axis.plot(mu_values,2.*cond_values_up,color='blue',marker='o')
#axis.plot(mu_extra,2.*cond_extra,color='red',marker='o')
#
##axis.set_xlim([-5,5])
#
#plt.show()

