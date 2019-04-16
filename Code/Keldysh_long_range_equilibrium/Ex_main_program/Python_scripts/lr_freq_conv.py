import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys
import One_particle_conductance_to_numpy as Opctn 

def load_conductance(s):
	liste = glob.glob(s)
	liste.sort()
	print(liste)
	cond_values_up = [] 
	cond_values_down = [] 
	mu_values =[]
	for s in liste:
		print(s)
		cond_values_up.append(np.real(ctn.load_as_np(s,"cond_up"))[0,:])
		cond_values_down.append(np.real(ctn.load_as_np(s,"cond_down"))[0,:])
		mu_values.append(np.real(ctn.load_as_np(s,"mu"))[0,:])
	
	mu_values = np.asarray(mu_values)
	cond_values_up = np.asarray(cond_values_up)
	cond_values_down = np.asarray(cond_values_down)
	return [mu_values, cond_values_up, cond_values_down]

[mu_val, cond_val_up, cond_val_down] = load_conductance('/naslx/projects/uh3o1/ri26yad/Ex_DATA_lrz/QPC_zero_tmp_lrz/Conductance/cond_ex_dsfRG_L_5_Lu_3_N_30_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_0_Vg_0.250000_h_0.000000_mu_*_T_0.000000_U0_0.500000_U1_0.300000_Xi_5.000000.mat')

[mu_val_Nfb500, cond_val_up_Nfb500, cond_val_down_Nfb500] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/Conductance/cond_X_L5_Lu3_N30_Nff1500_NfbP500_NfbX500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.30_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat')

[mu_val_Nfb250, cond_val_up_Nfb250, cond_val_down_Nfb250] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/Conductance/cond_X_L5_Lu3_N30_Nff1500_NfbP250_NfbX250_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.30_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat')

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(mu_val,cond_val_up+cond_val_down,color='blue',marker='o')
axis.plot(mu_val_Nfb250,cond_val_up_Nfb250+cond_val_down_Nfb250,color='green',marker='o')
axis.plot(mu_val_Nfb500,cond_val_up_Nfb500+cond_val_down_Nfb500,color='red',marker='o')

#axis.set_xlim([-5,5])
plt.savefig('lr_freq_conv.pdf', format = 'pdf')

plt.show()



