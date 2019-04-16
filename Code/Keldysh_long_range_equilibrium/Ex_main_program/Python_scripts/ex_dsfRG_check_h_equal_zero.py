import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str




#Load normal Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_save_tests/Schar_ex_dsfRG/"
s_file="ex_dsfRG_L_1_Lu_1_N_1_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_2_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000_ap_0.000100_ax_0.000100_as_0.000100_tol_0.000001_Lini_100000.000000_Lfin_0.000000_nodes_4.mat"

s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
L = np.real(ctn.load_as_np(s,"L"))[0,0]
L = L.astype(int)
wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
Lp_structure = Lp_structure.astype(int)
Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
Lx_structure = Lx_structure.astype(int)

NfbP = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
NfbP = NfbP.astype(int)
NfbX = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
NfbX = NfbX.astype(int)

h = np.real(ctn.load_as_np(s,"h"))[0,0]

[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

#Load h_equal_zero Flow-Data:
#s_folder="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_save_tests_h_equal_zero/Schar_ex_dsfRG/"
#s_file="ex_dsfRG_L_1_Lu_1_N_1_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_2_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000_ap_0.000100_ax_0.000100_as_0.000100_tol_0.000001_Lini_100000.000000_Lfin_0.000000_nodes_4.mat"

s_folder="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_save_tests_h_equal_zero_modified/Schar_ex_dsfRG/"
s_file="ex_dsfRG_L_1_Lu_1_N_1_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_2_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000.mat"

#s_folder="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_save_tests_modified/Schar_ex_dsfRG/"
#s_file="ex_dsfRG_L_1_Lu_1_N_1_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_2_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000.mat"

s=s_folder + s_file

N_hz = np.real(ctn.load_as_np(s,"N"))[0,0]
N_hz = N.astype(int)
L_hz = np.real(ctn.load_as_np(s,"L"))[0,0]
L_hz = L.astype(int)
wbP_hz = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_hz = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf_hz = np.real(ctn.load_as_np(s,"wf"))[0,:]
Lp_structure_hz = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
Lp_structure_hz = Lp_structure.astype(int)
Lx_structure_hz = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
Lx_structure_hz = Lx_structure.astype(int)

NfbP_hz = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
NfbP_hz = NfbP.astype(int)
NfbX_hz = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
NfbX_hz = NfbX.astype(int)

h_hz = np.real(ctn.load_as_np(s,"h"))[0,0]

[ERetu_hz, ERetd_hz, aPuu_dyn_hz, aPuu_stat_hz, aPdd_dyn_hz, aPdd_stat_hz, aPud_dyn_hz, aPud_stat_hz, aXud_dyn_hz, aXud_stat_hz, aDuu_dyn_hz, aDuu_stat_hz, aDdd_dyn_hz, aDdd_stat_hz, aDud_dyn_hz, aDud_stat_hz] = Ex_Vertex.load(s,"gamma")

print("begin comparison")
print(np.amax(np.absolute( ERetu - ERetu_hz)))
print(np.amax(np.absolute( ERetd - ERetd_hz)))
print(Ex_freq_str.abs_dyn(aPuu_dyn - aPuu_dyn_hz, Lp_structure))
print(Ex_freq_str.abs_stat(aPuu_stat - aPuu_stat_hz, L))
print(Ex_freq_str.abs_dyn(aPdd_dyn - aPdd_dyn_hz, Lp_structure))
print(Ex_freq_str.abs_stat(aPdd_stat - aPdd_stat_hz, L))
print(Ex_freq_str.abs_dyn(aPud_dyn - aPud_dyn_hz, Lp_structure))
print(Ex_freq_str.abs_stat(aPud_stat - aPud_stat_hz, L))
print(Ex_freq_str.abs_dyn(aXud_dyn - aXud_dyn_hz, Lx_structure))
print(Ex_freq_str.abs_stat(aXud_stat - aXud_stat_hz, L))
print(Ex_freq_str.abs_dyn(aDuu_dyn - aDuu_dyn_hz, Lx_structure))
print(Ex_freq_str.abs_stat(aDuu_stat - aDuu_stat_hz, L))
print(Ex_freq_str.abs_dyn(aDdd_dyn - aDdd_dyn_hz, Lx_structure))
print(Ex_freq_str.abs_stat(aDdd_stat - aDdd_stat_hz, L))
print(Ex_freq_str.abs_dyn(aDud_dyn - aDud_dyn_hz, Lx_structure))
print(Ex_freq_str.abs_stat(aDud_stat - aDud_stat_hz, L))
print("end")

#print(np.amax(np.absolute( aPuu_stat - aPuu_stat_hz)))
#print(np.amax(np.absolute( aPdd_dyn - aPdd_dyn_hz)))
#print(np.amax(np.absolute( aPdd_stat - aPdd_stat_hz)))
#print(np.amax(np.absolute( aPud_dyn - aPud_dyn_hz)))
#print(np.amax(np.absolute( aPud_stat - aPud_stat_hz)))
#print(np.amax(np.absolute( aXud_dyn - aXud_dyn_hz)))
#print(np.amax(np.absolute( aXud_stat - aXud_stat_hz)))
#print(np.amax(np.absolute( aDuu_dyn - aDuu_dyn_hz)))
#print(np.amax(np.absolute( aDuu_stat - aDuu_stat_hz)))
#print(np.amax(np.absolute( aDdd_dyn - aDdd_dyn_hz)))
#print(np.amax(np.absolute( aDdd_stat - aDdd_stat_hz)))
#print(np.amax(np.absolute( aDud_dyn - aDud_dyn_hz)))
#print(np.amax(np.absolute( aDud_stat - aDud_stat_hz)))



#fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
##	#axis.plot(wbX_old,np.real(aDuu_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbX_old,np.imag(aDuu_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='orange',linestyle='--',marker='o')
##	
##	#axis.plot(wbP_old,np.real(aPud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbP_old,np.imag(aPud_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbP,np.real(aPud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbP,np.imag(aPud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
##	
##	#axis.plot(wbX_old,np.real(aXud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbX_old,np.imag(aXud_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbX,np.real(aXud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
##	
#axis.plot(wf_old,np.real(ERetu_old[0,:,N,N]),color='blue',marker='o')
#axis.plot(wf_old,np.imag(ERetu_old[0,:,N,N]),color='red',marker='o')
#axis.plot(wf,np.real(ERetu[0,:,N,N]),color='black',linestyle='--',marker='o')
#axis.plot(wf,np.imag(ERetu[0,:,N,N]),color='orange',linestyle='--',marker='o')
#
#axis.set_xlim([-5,5])
#
#plt.show()


