import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

foldername="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test/N4_QPC/OMP_standard/"
file_s2 ="dsfRG_mpi_with_paid_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order2_number_of_nodes_5.mat" 
file_s4 ="dsfRG_mpi_with_paid_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order4_number_of_nodes_5.mat" 
file_s8 ="dsfRG_mpi_with_paid_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order8_number_of_nodes_5.mat" 
file_s16 ="dsfRG_mpi_with_paid_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order16_number_of_nodes_5.mat" 
s2 = foldername + file_s2
s4 = foldername + file_s4
s8 = foldername + file_s8
s16 = foldername + file_s16

#Check accuracy:
gamma_data_short_str_s2 = ctn.load_as_np(s2,"gamma_data_short_str")
gamma_data_short_str_s4 = ctn.load_as_np(s4,"gamma_data_short_str")
gamma_data_short_str_s8 = ctn.load_as_np(s4,"gamma_data_short_str")
err_short_str= np.amax(abs(gamma_data_short_str_s2 - gamma_data_short_str_s4))
print(err_short_str)
Puu_s2 = ctn.load_as_np(s2,"gamma_data_Puu")
Pdd_s2 = ctn.load_as_np(s2,"gamma_data_Pdd")
Pud_s2 = ctn.load_as_np(s2,"gamma_data_Pud")
Xud_s2 = ctn.load_as_np(s2,"gamma_data_Xud")
Duu_s2 = ctn.load_as_np(s2,"gamma_data_Duu")
Ddd_s2 = ctn.load_as_np(s2,"gamma_data_Ddd")
Dud_s2 = ctn.load_as_np(s2,"gamma_data_Dud")

Puu_s4= ctn.load_as_np(s4,"gamma_data_Puu")
Pdd_s4= ctn.load_as_np(s4,"gamma_data_Pdd")
Pud_s4= ctn.load_as_np(s4,"gamma_data_Pud")
Xud_s4= ctn.load_as_np(s4,"gamma_data_Xud")
Duu_s4= ctn.load_as_np(s4,"gamma_data_Duu")
Ddd_s4= ctn.load_as_np(s4,"gamma_data_Ddd")
Dud_s4= ctn.load_as_np(s4,"gamma_data_Dud")

Puu_diff = np.amax(abs(Puu_s4-Puu_s2))
Pdd_diff = np.amax(abs(Pdd_s4-Pdd_s2))
Pud_diff = np.amax(abs(Pud_s4-Pud_s2))
Xud_diff = np.amax(abs(Xud_s4-Xud_s2)) 
Duu_diff = np.amax(abs(Duu_s4-Duu_s2)) 
Ddd_diff = np.amax(abs(Ddd_s4-Ddd_s2)) 
Dud_diff = np.amax(abs(Dud_s4-Dud_s2)) 
all_diff=[Puu_diff, Pdd_diff, Pud_diff, Xud_diff, Duu_diff, Ddd_diff, Dud_diff]
print(all_diff)

error_total = max(np.amax(all_diff),err_short_str)
print(error_total)



#Check performance:
performance_track_s2 = ctn.load_as_np(s2,"performance_track");
comp_time_s2 = ctn.load_as_np(s2,"computation_time");
print(comp_time_s2)

performance_track_s4 = ctn.load_as_np(s4,"performance_track");
comp_time_s4 = ctn.load_as_np(s4,"computation_time");
print(comp_time_s4)

performance_track_s8 = ctn.load_as_np(s8,"performance_track");
comp_time_s8 = ctn.load_as_np(s8,"computation_time");
print(comp_time_s8)

performance_track_s16 = ctn.load_as_np(s16,"performance_track");
comp_time_s16 = ctn.load_as_np(s16,"computation_time");
print(comp_time_s16)


fig, axis = plt.subplots(nrows = 1, ncols = 3, figsize = (15,10))
axis[0].loglog(performance_track_s2[0,:,0,0], abs(performance_track_s2[0,:,0,0] - performance_track_s4[0,:,0,0]),marker='o')
axis[1].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,1],marker='o',color='blue')
#axis[1].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,4],marker='o',color='green')
#axis[1].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,7],marker='o',color='red')
#axis[1].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,1],marker='x',color='blue',linestyle='--')
#axis[1].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,4],marker='x',color='green',linestyle='--')
#axis[1].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,7],marker='x',color='red',linestyle='--')
axis[1].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,1],marker='D',color='pink',linestyle='--')
#axis[1].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,4],marker='D',color='green',linestyle='--')
#axis[1].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,7],marker='D',color='red',linestyle='--')

axis[2].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,2],marker='o',color='blue')
#axis[2].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,5],marker='o',color='green')
#axis[2].semilogx(performance_track_s2[0,:,0,0], performance_track_s2[0,:,0,8],marker='o',color='red')
#axis[2].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,2],marker='x',color='blue',linestyle='--')
#axis[2].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,5],marker='x',color='green',linestyle='--')
#axis[2].semilogx(performance_track_s4[0,:,0,0], performance_track_s4[0,:,0,8],marker='x',color='red',linestyle='--')
axis[2].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,2],marker='D',color='pink',linestyle='--')
#axis[2].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,5],marker='D',color='green',linestyle='--')
#axis[2].semilogx(performance_track_s8[0,:,0,0], performance_track_s8[0,:,0,8],marker='D',color='red',linestyle='--')
plt.show()
