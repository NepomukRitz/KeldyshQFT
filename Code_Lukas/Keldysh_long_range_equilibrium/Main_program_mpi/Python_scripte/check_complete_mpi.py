import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

#s1="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi/dsfRG_mpi_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"
#s2="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi/dsfRG_mpi_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_2.mat"
	
s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test/N15_QPC/Standard_code/"
s1_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (1)
s2_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (4)
s1=s_folder + s1_file
s2=s_folder + s2_file

gamma_data_short_str_complete_mpi = ctn.load_as_np(s1,"gamma_data_short_str")
gamma_data_short_str = ctn.load_as_np(s2,"gamma_data_short_str")
wf= ctn.load_as_np(s1,"wf");
wf=wf[0,:]

err_short_str= np.amax(abs(gamma_data_short_str_complete_mpi - gamma_data_short_str))
print(err_short_str)

Puu_complete_mpi = ctn.load_as_np(s1,"gamma_data_Puu")
Pdd_complete_mpi = ctn.load_as_np(s1,"gamma_data_Pdd")
Pud_complete_mpi = ctn.load_as_np(s1,"gamma_data_Pud")
Xud_complete_mpi = ctn.load_as_np(s1,"gamma_data_Xud")
Duu_complete_mpi = ctn.load_as_np(s1,"gamma_data_Duu")
Ddd_complete_mpi = ctn.load_as_np(s1,"gamma_data_Ddd")
Dud_complete_mpi = ctn.load_as_np(s1,"gamma_data_Dud")

Puu = ctn.load_as_np(s2,"gamma_data_Puu")
Pdd = ctn.load_as_np(s2,"gamma_data_Pdd")
Pud = ctn.load_as_np(s2,"gamma_data_Pud")
Xud = ctn.load_as_np(s2,"gamma_data_Xud")
Duu = ctn.load_as_np(s2,"gamma_data_Duu")
Ddd = ctn.load_as_np(s2,"gamma_data_Ddd")
Dud = ctn.load_as_np(s2,"gamma_data_Dud")

Puu_diff = np.amax(abs(Puu-Puu_complete_mpi))
Pdd_diff = np.amax(abs(Pdd-Pdd_complete_mpi))
Pud_diff = np.amax(abs(Pud-Pud_complete_mpi))
Xud_diff = np.amax(abs(Xud-Xud_complete_mpi)) 
Duu_diff = np.amax(abs(Duu-Duu_complete_mpi)) 
Ddd_diff = np.amax(abs(Ddd-Ddd_complete_mpi)) 
Dud_diff = np.amax(abs(Dud-Dud_complete_mpi)) 

all_diff=[Puu_diff, Pdd_diff, Pud_diff, Xud_diff, Duu_diff, Ddd_diff, Dud_diff]
print(all_diff)
#
#
fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
#axis.plot(wf,np.imag(gamma_data_short_str_complete_mpi[0,3,0,:,4,4]))
#axis.plot(wf,np.imag(gamma_data_short_str[0,3,0,:,4,4]))
axis.plot(np.imag(gamma_data_short_str[0,7,0,:,4,4] - gamma_data_short_str_complete_mpi[0,7,0,:,4,4]),marker='o')
print(np.imag(gamma_data_short_str[0,7,0,-1,4,4] - gamma_data_short_str_complete_mpi[0,7,0,-1,4,4]))
#axis.set_xlim([-5,5])
#axis[0].imshow(np.real(Xud_complete_mpi))
#axis[1].imshow(np.real(Xud))
#axis.imshow(np.real(Puu), vmin=-1e-5, vmax=1e-5)
plt.show()
