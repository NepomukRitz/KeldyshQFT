import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import glob 
import sys


#Old mpi-data:

s="/p/scratch/chmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0/Schar_dsfRG_mpi_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/dsfRG_mpi_L0_Lu0_N30_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.000000_U0_0.650000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_8.mat"

wf_old = np.real(ctn.load_as_np(s,"wf"))[0,:]
gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
print(np.shape(gamma_data_short_str))
self_old=gamma_data_short_str[0,7,0,:,:,:]

#New data:

s="/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp/ex_dsfRG_L_0_Lu_0_N_30_Nff_1500_NfbP_1500_NfbX_1500_npre_30000_NL_0_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.000000_U0_0.650000_U1_0.000000_Xi_5.000000.mat"

wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

print(np.amax(np.absolute(wf_old-wf)))

diff_self=ERetu - self_old
print("diff_self=",np.amax(np.absolute(diff_self)))

diff=np.empty([np.size(wf)])
for i in range(np.size(wf)):
	diff[i] = np.amax(np.absolute(ERetu[0,i,:,:] - self_old[i,:,:]))
	

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(wf_old,np.imag(self_old[:,30,30]),color='blue',marker='o')
axis.plot(wf,np.imag(ERetu[0,:,30,30]),color='red',marker='o')
axis.plot(wf_old,np.real(self_old[:,30,30]),color='green',marker='o')
axis.plot(wf,np.real(ERetu[0,:,30,30]),color='orange',marker='o')
#axis.plot(wf,diff,color='orange',marker='o')

axis.set_xlim([-5,5])

plt.show()

