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


[mu_val_Nfb1500, cond_val_up_Nfb1500, cond_val_down_Nfb1500] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
[mu_val_Nfb500, cond_val_up_Nfb500, cond_val_down_Nfb500] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP500_NfbX500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
[mu_val_Nfb250, cond_val_up_Nfb250, cond_val_down_Nfb250] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
[mu_val_Nfb125, cond_val_up_Nfb125, cond_val_down_Nfb125] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP125_NfbX125_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')

s="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/X_L0_Lu0_N30_Nff1500_NfbP125_NfbX125_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4500_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
mu_list = np.linspace(-1.6,-1.4,50)
cond_non_int = Opctn.cond_combined(s,mu_list)


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(mu_val_Nfb1500,cond_val_up_Nfb1500+cond_val_down_Nfb1500,color='blue',marker='o')
axis.plot(mu_val_Nfb500,cond_val_up_Nfb500+cond_val_down_Nfb500,color='red',marker='o')
axis.plot(mu_val_Nfb250,cond_val_up_Nfb250+cond_val_down_Nfb250,color='orange',marker='o')
axis.plot(mu_val_Nfb125,cond_val_up_Nfb125+cond_val_down_Nfb125,color='magenta',marker='o')
axis.plot(mu_list,cond_non_int,color='black',marker='o')

plt.savefig('sr_freq_conv.pdf', format = 'pdf')

plt.show()



