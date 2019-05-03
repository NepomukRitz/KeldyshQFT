import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

s="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi/dsfRG_mpi_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_2.mat"
s_with_paid="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi_with_paid/dsfRG_mpi_with_paid_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_2.mat"

#Check accuracy:
gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
gamma_data_short_str_with_paid = ctn.load_as_np(s_with_paid,"gamma_data_short_str")
Nges = ctn.load_as_np(s,"Nges")[0,0]
size_syma=(Nges*Nges + Nges)/2
print(size_syma)
err_short_str= np.amax(abs(gamma_data_short_str - gamma_data_short_str_with_paid))
print(err_short_str)
Puu = ctn.load_as_np(s,"gamma_data_Puu")
Pdd = ctn.load_as_np(s,"gamma_data_Pdd")
Pud = ctn.load_as_np(s,"gamma_data_Pud")
Xud = ctn.load_as_np(s,"gamma_data_Xud")
Duu = ctn.load_as_np(s,"gamma_data_Duu")
Ddd = ctn.load_as_np(s,"gamma_data_Ddd")
Dud = ctn.load_as_np(s,"gamma_data_Dud")

Puu_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Puu")
Pdd_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Pdd")
Pud_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Pud")
Xud_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Xud")
Duu_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Duu")
Ddd_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Ddd")
Dud_with_paid= ctn.load_as_np(s_with_paid,"gamma_data_Dud")

Puu_diff = np.amax(abs(Puu_with_paid-Puu))
Pdd_diff = np.amax(abs(Pdd_with_paid-Pdd))
Pud_diff = np.amax(abs(Pud_with_paid-Pud))
Xud_diff = np.amax(abs(Xud_with_paid-Xud)) 
Duu_diff = np.amax(abs(Duu_with_paid-Duu)) 
Ddd_diff = np.amax(abs(Ddd_with_paid-Ddd)) 
Dud_diff = np.amax(abs(Dud_with_paid-Dud)) 
all_diff=[Puu_diff, Pdd_diff, Pud_diff, Xud_diff, Duu_diff, Ddd_diff, Dud_diff]
print(all_diff)

error_total = max(np.amax(all_diff),err_short_str)
print(error_total)



#Check performance:
performance_track = ctn.load_as_np(s,"performance_track");
comp_time = ctn.load_as_np(s,"computation_time");
print(comp_time)

performance_track_with_paid = ctn.load_as_np(s_with_paid,"performance_track");
comp_time_with_paid = ctn.load_as_np(s_with_paid,"computation_time");
print(comp_time_with_paid)


fig, axis = plt.subplots(nrows = 1, ncols = 3, figsize = (15,10))
axis[0].loglog(performance_track[0,:,0,0], abs(performance_track[0,:,0,0] - performance_track_with_paid[0,:,0,0]),marker='o')
axis[1].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,1],marker='o',color='blue')
axis[1].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,1],marker='o',color='blue',linestyle='--')
axis[1].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,4],marker='o',color='green')
axis[1].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,4],marker='o',color='green',linestyle='--')
axis[1].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,7],marker='o',color='red')
axis[1].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,7],marker='o',color='red',linestyle='--')

axis[2].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,2]*(Nges*Nges + Nges)/2,marker='o',color='blue')
axis[2].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,2],marker='o',color='blue',linestyle='--')
axis[2].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,5]*(Nges*Nges + Nges)/2,marker='o',color='green')
axis[2].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,5],marker='o',color='green',linestyle='--')
axis[2].semilogx(performance_track[0,:,0,0], performance_track[0,:,0,8]*(Nges*Nges + Nges)/2,marker='o',color='red')
axis[2].semilogx(performance_track_with_paid[0,:,0,0], performance_track_with_paid[0,:,0,8],marker='o',color='red',linestyle='--')
plt.show()
