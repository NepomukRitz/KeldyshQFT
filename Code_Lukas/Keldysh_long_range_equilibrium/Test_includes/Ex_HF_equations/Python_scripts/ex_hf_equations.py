import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex




#Load Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG/"
s_file="ex_dsfRG_add_katanin_L0_Lu0_N2_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_0.000000_T_0.000000_U0_0.100000_U1_0.000000_Xi_5.000000_NL_full_0_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"

s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

#Load HF-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_HF_equations/"
s_file="ex_hf_equations_L0_Lu0_N2_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_0.000000_T_0.000000_U0_0.100000_U1_0.000000_Xi_5.000000_NL_full_0_Lambda_ini_0.000000_Lambda_fin_0.000000_number_of_nodes_1.mat"

s=s_folder + s_file

wbP_hf = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_hf = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf_hf = np.real(ctn.load_as_np(s,"wf"))[0,:]
[ERetu_hf, ERetd_hf, aPuu_dyn_hf, aPuu_stat_hf, aPdd_dyn_hf, aPdd_stat_hf, aPud_dyn_hf, aPud_stat_hf, aXud_dyn_hf, aXud_stat_hf, aDuu_dyn_hf, aDuu_stat_hf, aDdd_dyn_hf, aDdd_stat_hf, aDud_dyn_hf, aDud_stat_hf] = Ex_Vertex.load(s,"gamma")

print("diff_ERetu=",np.amax(np.absolute(ERetu_hf - ERetu)))
print("diff_ERetd=",np.amax(np.absolute(ERetd_hf - ERetd)))


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#	#axis.plot(wbX_hf,np.real(aDuu_dyn_hf[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_hf,np.imag(aDuu_dyn_hf[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbP_hf,np.real(aPud_dyn_hf[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbP_hf,np.imag(aPud_dyn_hf[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbP,np.real(aPud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbP,np.imag(aPud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbX_hf,np.real(aXud_dyn_hf[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_hf,np.imag(aXud_dyn_hf[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aXud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	
axis.plot(wf_hf,np.real(ERetu_hf[0,:,N,N]),color='blue',marker='o')
axis.plot(wf_hf,np.imag(ERetu_hf[0,:,N,N]),color='red',marker='o')
axis.plot(wf,np.real(ERetu[0,:,N,N]),color='black',linestyle='--',marker='o')
axis.plot(wf,np.imag(ERetu[0,:,N,N]),color='orange',linestyle='--',marker='o')

axis.set_xlim([-5,5])

plt.show()


