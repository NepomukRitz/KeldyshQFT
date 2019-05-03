import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys

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

#[mu_val, cond_val_up, cond_val_down] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
##
#[mu_val_NL5, cond_val_up_NL5, cond_val_down_NL5] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL5_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no10.mat')
##
#[mu_val_NL2, cond_val_up_NL2, cond_val_down_NL2] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL2_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat')
##
#[mu_val_NL10, cond_val_up_NL10, cond_val_down_NL10] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL10_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no8.mat')
#[mu_val_NL10nle, cond_val_up_NL10nle, cond_val_down_NL10nle] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL10_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no8.mat')




[mu_val_Nfb1500, cond_val_up_Nfb1500, cond_val_down_Nfb1500] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
[mu_val_Nfb500, cond_val_up_Nfb500, cond_val_down_Nfb500] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP500_NfbX500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')
[mu_val_Nfb250, cond_val_up_Nfb250, cond_val_down_Nfb250] = load_conductance('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/Conductance/cond_X_L0_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#axis.plot(mu_val,cond_val_up+cond_val_down,color='blue',marker='o')
#axis.plot(mu_val_NL10nle,cond_val_up_NL10nle+cond_val_down_NL10nle,color='red',marker='o')
#axis.plot(mu_val_NL2,cond_val_up_NL2+cond_val_down_NL2,color='green',marker='o')
#axis.plot(mu_val_NL5,cond_val_up_NL5+cond_val_down_NL5,color='red',marker='o')
#axis.plot(mu_val_NL10,cond_val_up_NL10+cond_val_down_NL10,color='orange',marker='o')

axis.plot(mu_val_Nfb1500,cond_val_up_Nfb1500+cond_val_down_Nfb1500,color='blue',marker='o')
axis.plot(mu_val_Nfb500,cond_val_up_Nfb500+cond_val_down_Nfb500,color='red',marker='o')
axis.plot(mu_val_Nfb250,cond_val_up_Nfb250+cond_val_down_Nfb250,color='orange',marker='o')

#axis.set_xlim([-5,5])

plt.show()


