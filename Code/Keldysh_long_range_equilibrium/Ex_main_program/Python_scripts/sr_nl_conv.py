import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys
from math import *

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

[mu_val_NL0, cond_val_up_NL0, cond_val_down_NL0] = load_conductance('/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp_no_lr_extrapolation/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat')

[mu_val_NL5, cond_val_up_NL5, cond_val_down_NL5] = load_conductance('/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp_no_lr_extrapolation/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL5_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat')

[mu_val_NL10, cond_val_up_NL10, cond_val_down_NL10] = load_conductance('/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp_no_lr_extrapolation/Conductance/cond_X_L5_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL10_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no8.mat')


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(mu_val_NL0,cond_val_up_NL0+cond_val_down_NL0,color='blue',marker='o')
axis.plot(mu_val_NL5,cond_val_up_NL5+cond_val_down_NL5,color='red',marker='o')
axis.plot(mu_val_NL10,cond_val_up_NL10+cond_val_down_NL10,color='orange',marker='o')

axis.set_xlim([-1.6,-1.35])
axis.set_ylim([0,1])
plt.savefig('sr_nl_conv.pdf', format = 'pdf')

plt.show()


