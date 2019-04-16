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





[mu_val_T0, cond_val_up_T0, cond_val_down_T0] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')

[mu_val_T0005, cond_val_up_T0005, cond_val_down_T0005] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.005000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')

#[mu_val_T0005_old, cond_val_up_T0005_old, cond_val_down_T0005_old] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation_old/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.005000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')



fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(mu_val_T0,cond_val_up_T0+cond_val_down_T0,color='blue',marker='o')
axis.plot(mu_val_T0005,0.5*cond_val_up_T0005+0.5*cond_val_down_T0005,color='red',marker='o')
#axis.plot(mu_val_T0005_old,0.5*cond_val_up_T0005_old+0.5*cond_val_down_T0005_old,color='red',linestyle='--')


#axis.set_xlim([-5,5])
plt.savefig('sr_finite_temp.pdf', format = 'pdf')

plt.show()


