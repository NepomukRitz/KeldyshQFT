import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

paid_order = 4
number_of_nodes = 8 

#s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test/N15_QPC/OMP_standard/"
s1_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/OMP_standard/"
s2_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/Standard_code/"

s1_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_%d.mat" % (paid_order,number_of_nodes)
s2_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (number_of_nodes)
s1=s1_folder + s1_file
s2=s2_folder + s2_file

performance_track_paid = ctn.load_as_np(s1,"performance_track");
computation_time_paid = ctn.load_as_np(s1,"computation_time");
computation_time_standard = ctn.load_as_np(s2,"computation_time");
performance_track_standard = ctn.load_as_np(s2,"performance_track");

print(performance_track_paid[0,:,0,0] - performance_track_standard[0,:,0,0])
print(computation_time_paid)
print(computation_time_standard)
print(computation_time_paid/computation_time_standard)

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))

axis[0].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,1], color='blue',marker='o')
axis[0].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,4], color='green',marker='o')
axis[0].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,7], color='red',marker='o')

axis[0].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,1], color='blue',marker='D')
axis[0].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,4], color='green',marker='D')
axis[0].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,7], color='red',marker='D')

axis[1].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,2]*31*31, color='blue',marker='o')
axis[1].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,5]*31*31, color='green',marker='o')
axis[1].loglog(performance_track_standard[0,:,0,0],performance_track_standard[0,:,0,8]*31*31, color='red',marker='o')

axis[1].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,2], color='blue',marker='D')
axis[1].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,5], color='green',marker='D')
axis[1].loglog(performance_track_paid[0,:,0,0],performance_track_paid[0,:,0,8], color='red',marker='D')

axis[0].set_xlabel(r'Flow parameter $\Lambda$')
axis[0].set_ylabel(r'Computation_time/seconds')

axis[0].text(1e-8,2.5,r'$\Sigma$')
axis[0].text(1e-8,1.3,r'$X$')
axis[0].text(1e-8,0.8,r'$P$')
axis[0].text(1e-6,3,r'Old integrator')
axis[0].text(1e-6,40,r'PAID')

axis[1].set_xlabel(r'Flow parameter $\Lambda$')
axis[1].set_ylabel(r'Integrand evaluations')
axis[1].text(1e-8,9e8,r'Old integrator')
axis[1].text(1e-8,2.5e7,r'PAID')
plt.savefig('qpc_performance_flow_resolved_8_nodes.pdf', format = 'pdf')
fig.suptitle(r'QPC, $31$ sites, $8$ Nodes')

plt.show()
