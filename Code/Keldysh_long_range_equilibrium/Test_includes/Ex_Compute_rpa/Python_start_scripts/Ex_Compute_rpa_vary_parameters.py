import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys

def imaginary_part(A,Nff,N,L_structure):
	value=np.empty([Nff])
	for i in range(Nff):
	 	L_inner = L_structure[i]
		value[i] = np.imag(A[i][L_inner,L_inner][N,N])
	return value

def real_part(A,Nff,N,L_structure):
	value=np.empty([Nff])
	for i in range(Nff):
	 	L_inner = L_structure[i]
		value[i] = np.real(A[i][L_inner,L_inner][N,N])
	return value

def absolute(A,Nff,N,L_structure):
	value=np.empty([Nff])
	for i in range(Nff):
	 	L_inner = L_structure[i]
		value[i] = np.absolute(A[i][L_inner,L_inner][N,N])
	return value

def load_vertex_schar(s,vertex_variable):
	liste = glob.glob(s)
	liste.sort()
	print(liste)
	N_values =[]
	NfbP_values =[]
	NfbX_values =[]
	NL_full_values =[]
	wbP_values =[]
	wbX_values =[]
	Lp_structure_values =[]
	Lx_structure_values =[]
	aPud_dyn_values = [] 
	aXud_dyn_values = [] 
	for s in liste:
		print(s)
		N = np.real(ctn.load_as_np(s,"N"))[0,0]
		N = N.astype(int)
		N_values.append(N)
		NfbP = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
		NfbP = NfbP.astype(int)
		NfbP_values.append(NfbP)
		NfbX = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
		NfbX = NfbX.astype(int)
		NfbX_values.append(NfbX)
		NL_full = np.real(ctn.load_as_np(s,"NL_full"))[0,0]
		NL_full = NL_full.astype(int)
		NL_full_values.append(NL_full)
		wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
		wbP_values.append(wbP)
		wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
		wbX_values.append(wbX)
		Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
		Lp_structure_values.append(Lp_structure)
		Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
		Lx_structure_values.append(Lx_structure)
		[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,vertex_variable)
		aPud_dyn_values.append(aPud_dyn)
		aXud_dyn_values.append(aXud_dyn)
	
	N_values = np.asarray(N_values)
	NfbP_values = np.asarray(NfbP_values)
	NfbX_values = np.asarray(NfbX_values)
	NL_full_values = np.asarray(NL_full_values)
	wbP_values = np.asarray(wbP_values)
	wbX_values = np.asarray(wbX_values)
	Lp_structure_values = np.asarray(Lp_structure_values)
	Lx_structure_values = np.asarray(Lx_structure_values)
	#aPud_dyn_values = np.asarray(aPud_dyn_values)
	#aXud_dyn_values = np.asarray(aXud_dyn_values)
	return [N_values, NfbP_values, NfbX_values, NL_full_values, wbP_values, wbX_values, Lp_structure_values, Lx_structure_values, aPud_dyn_values, aXud_dyn_values]
	 	

#Load RPA-Data:	
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Compute_rpa/"
s_file="RPA_L2_Lu0_N5_Nff1500_NfbP1500_NfbX1500_pre30000_NL*_Vg0.2500_h0.000000_mu-1.4750_T0.000000_Uc0.30_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04.mat"
vertex_variable = "gamma_rpa"


s=s_folder + s_file
[N_values, NfbP_values, NfbX_values, NL_full_values, wbP_values, wbX_values, Lp_structure_values, Lx_structure_values, aPud_dyn_values, aXud_dyn_values] = load_vertex_schar(s,vertex_variable)

#wbP_rpa = np.real(ctn.load_as_np(s,"wbP"))[0,:]
#wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
#[ERetu_rpa, ERetd_rpa, aPuu_dyn_rpa, aPuu_stat_rpa, aPdd_dyn_rpa, aPdd_stat_rpa, aPud_dyn_rpa, aPud_stat_rpa, aXud_dyn_rpa, aXud_stat_rpa, aDuu_dyn_rpa, aDuu_stat_rpa, aDdd_dyn_rpa, aDdd_stat_rpa, aDud_dyn_rpa, aDud_stat_rpa] = Ex_Vertex.load(s,"gamma_rpa")
#
#N_rpa = np.real(ctn.load_as_np(s,"N"))[0,0]
#N_rpa = N_rpa.astype(int)
#L_rpa = np.real(ctn.load_as_np(s,"L"))[0,0]
#L_rpa = L_rpa.astype(int)
#wbP_rpa = np.real(ctn.load_as_np(s,"wbP"))[0,:]
#wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
#wf_rpa = np.real(ctn.load_as_np(s,"wf"))[0,:]
#Lp_structure_rpa = Lp_structure_rpa.astype(int)
#Lx_structure_rpa = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
#Lx_structure_rpa = Lx_structure_rpa.astype(int)
#NfbX_rpa = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
#NfbX_rpa = NfbX_rpa.astype(int)
#
#
#
fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (6,3))

axis[0].plot(wbP_values[0],real_part(aPud_dyn_values[0],NfbP_values[0],N_values[0],Lp_structure_values[0]),color='blue',marker='o')
axis[0].plot(wbP_values[1],real_part(aPud_dyn_values[1],NfbP_values[1],N_values[1],Lp_structure_values[1]),color='red',marker='o')

axis[0].set_xlim([-6,6])

plt.show()



