import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex




#Load Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG/"
s_file="ex_dsfRG_L0_Lu0_N2_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.020000_mu_0.000000_T_0.010000_U0_0.100000_U1_0.000000_Xi_5.000000_NL_full_0_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"

s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

#Load Old-Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG/"
s_file="ex_dsfRG_vgl_flow_old_L0_Lu0_N2_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.020000_mu_0.000000_T_0.010000_U0_0.100000_U1_0.000000_Xi_5.000000_NL_full_0_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"

s=s_folder + s_file

wbP_old = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_old = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf_old = np.real(ctn.load_as_np(s,"wf"))[0,:]
[ERetu_old, ERetd_old, aPuu_dyn_old, aPuu_stat_old, aPdd_dyn_old, aPdd_stat_old, aPud_dyn_old, aPud_stat_old, aXud_dyn_old, aXud_stat_old, aDuu_dyn_old, aDuu_stat_old, aDdd_dyn_old, aDdd_stat_old, aDud_dyn_old, aDud_stat_old] = Ex_Vertex.load(s,"gamma")

print("diff_ERetu=",np.amax(np.absolute(ERetu_old - ERetu)))
print("diff_ERetd=",np.amax(np.absolute(ERetd_old - ERetd)))
print("diff_pud=",np.amax(np.absolute(aPud_dyn_old - aPud_dyn)))
print("diff_xud=",np.amax(np.absolute(aXud_dyn_old - aXud_dyn)))
print("diff_duu=",np.amax(np.absolute(aDuu_dyn_old - aDuu_dyn)))
print("diff_ddd=",np.amax(np.absolute(aDdd_dyn_old - aDdd_dyn)))


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#	#axis.plot(wbX_old,np.real(aDuu_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_old,np.imag(aDuu_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbP_old,np.real(aPud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbP_old,np.imag(aPud_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbP,np.real(aPud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbP,np.imag(aPud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbX_old,np.real(aXud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_old,np.imag(aXud_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aXud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	
axis.plot(wf_old,np.real(ERetu_old[0,:,N,N]),color='blue',marker='o')
axis.plot(wf_old,np.imag(ERetu_old[0,:,N,N]),color='red',marker='o')
axis.plot(wf,np.real(ERetu[0,:,N,N]),color='black',linestyle='--',marker='o')
axis.plot(wf,np.imag(ERetu[0,:,N,N]),color='orange',linestyle='--',marker='o')

axis.set_xlim([-5,5])

plt.show()


