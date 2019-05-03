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

def load_vertex_maxima(s):
	liste = glob.glob(s)
	liste.sort()
	print(liste)
	aPuu_dyn_max = [] 
	aPuu_stat_max = [] 
	aPdd_dyn_max = [] 
	aPdd_stat_max = [] 
	aPud_dyn_max = [] 
	aPud_stat_max = [] 
	aXud_dyn_max = [] 
	aXud_stat_max = [] 
	aDuu_dyn_max = [] 
	aDuu_stat_max = [] 
	aDdd_dyn_max = [] 
	aDdd_stat_max = [] 
	aDud_dyn_max = [] 
	aDud_stat_max = [] 
	mu_values =[]
	for s in liste:
		print(s)
		[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")
		mu_values.append(np.real(ctn.load_as_np(s,"mu"))[0,:])
		L = np.real(ctn.load_as_np(s,"L"))[0,0]
		L = L.astype(int)
		Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
		Lp_structure = Lp_structure.astype(int)
		Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
		Lx_structure = Lx_structure.astype(int)
		print(Lp_structure[0])
		print(Lp_structure[0]+77)
		aPuu_dyn_max.append(Ex_freq_str.abs_dyn(aPuu_dyn,Lp_structure))
		aPuu_stat_max.append(Ex_freq_str.abs_stat(aPuu_stat,L))
		aPdd_dyn_max.append(Ex_freq_str.abs_dyn(aPdd_dyn,Lp_structure))
		aPdd_stat_max.append(Ex_freq_str.abs_stat(aPdd_stat,L))
		aPud_dyn_max.append(Ex_freq_str.abs_dyn(aPud_dyn,Lp_structure))
		aPud_stat_max.append(Ex_freq_str.abs_stat(aPud_stat,L))
		aXud_dyn_max.append(Ex_freq_str.abs_dyn(aXud_dyn,Lx_structure))
		aXud_stat_max.append(Ex_freq_str.abs_stat(aXud_stat,L))
		aDuu_dyn_max.append(Ex_freq_str.abs_dyn(aDuu_dyn,Lx_structure))
		aDuu_stat_max.append(Ex_freq_str.abs_stat(aDuu_stat,L))
		aDdd_dyn_max.append(Ex_freq_str.abs_dyn(aDdd_dyn,Lx_structure))
		aDdd_stat_max.append(Ex_freq_str.abs_stat(aDdd_stat,L))
		aDud_dyn_max.append(Ex_freq_str.abs_dyn(aDud_dyn,Lx_structure))
		aDud_stat_max.append(Ex_freq_str.abs_stat(aDud_stat,L))
	
	mu_values = np.asarray(mu_values)
	aPuu_dyn_max = np.asarray(aPuu_dyn_max)
	aPuu_stat_max = np.asarray(aPuu_stat_max)
	aPdd_dyn_max = np.asarray(aPdd_dyn_max)
	aPdd_stat_max = np.asarray(aPdd_stat_max)
	aPud_dyn_max = np.asarray(aPud_dyn_max)
	aPud_stat_max = np.asarray(aPud_stat_max)
	aXud_dyn_max = np.asarray(aXud_dyn_max)
	aXud_stat_max = np.asarray(aXud_stat_max)
	aDuu_dyn_max = np.asarray(aDuu_dyn_max)
	aDuu_stat_max = np.asarray(aDuu_stat_max)
	aDdd_dyn_max = np.asarray(aDdd_dyn_max)
	aDdd_stat_max = np.asarray(aDdd_stat_max)
	aDud_dyn_max = np.asarray(aDud_dyn_max)
	aDud_stat_max = np.asarray(aDud_stat_max)
	return [mu_values, aPuu_dyn_max, aPuu_stat_max, aPdd_dyn_max, aPdd_stat_max, aPud_dyn_max, aPud_stat_max, aXud_dyn_max, aXud_stat_max, aDuu_dyn_max, aDuu_stat_max, aDdd_dyn_max, aDdd_stat_max, aDud_dyn_max, aDud_stat_max]


[mu_values, aPuu_dyn_max, aPuu_stat_max, aPdd_dyn_max, aPdd_stat_max, aPud_dyn_max, aPud_stat_max, aXud_dyn_max, aXud_stat_max, aDuu_dyn_max, aDuu_stat_max, aDdd_dyn_max, aDdd_stat_max, aDud_dyn_max, aDud_stat_max] = load_vertex_maxima('/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu*_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat')


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(mu_values,aPud_dyn_max,color='blue',marker='o')
#axis.plot(mu_val_Lu3,cond_val_up_Lu3+cond_val_down_Lu3,color='blue',marker='o')
#
##axis.plot(mu_val,cond_val_up+cond_val_down,color='blue',marker='o')
##axis.plot(mu_val_NL2,cond_val_up_NL2+cond_val_down_NL2,color='green',marker='o')
##axis.plot(mu_val_NL5,cond_val_up_NL5+cond_val_down_NL5,color='red',marker='o')
##axis.plot(mu_val_NL10,cond_val_up_NL10+cond_val_down_NL10,color='orange',marker='o')
#
##axis.set_xlim([-5,5])
#
plt.show()


