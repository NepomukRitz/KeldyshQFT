import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

#s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/OMP_standard/"
s_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/OMP_standard/"

paid_order_list=[2,4,6,8,16]
#paid_order_list=[2,4]


#Preload to determine dimensions:
s_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_8.mat" % (paid_order_list[0])
s=s_folder + s_file
Nges = ctn.load_as_np(s,"Nges")[0,0].astype(int);
L = ctn.load_as_np(s,"L")[0,0];
wf = ctn.load_as_np(s,"wf")[0,:];
Puu = ctn.load_as_np(s,"gamma_data_Puu")
Block_dim = np.shape(Puu)[0]
Nff_real = len(wf)

gamma_data_short_str_list=np.empty([len(paid_order_list),1,9, 1, Nff_real,Nges,Nges],dtype=complex)
gamma_data_long_str_list=np.empty([len(paid_order_list),7,Block_dim,Block_dim])

comp_time_list =np.empty(len(paid_order_list))

#Load all data:
for paid_order_ind in range(len(paid_order_list)):
	paid_order = paid_order_list[paid_order_ind]
	s_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_8.mat" % (paid_order)
	s=s_folder + s_file
	gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
	gamma_data_short_str_list[paid_order_ind,:,:,:,:,:,:] = gamma_data_short_str[:,:,:,:,:,:]
	Puu = ctn.load_as_np(s,"gamma_data_Puu")
	Pdd = ctn.load_as_np(s,"gamma_data_Pdd")
	Pud = ctn.load_as_np(s,"gamma_data_Pud")
	Xud = ctn.load_as_np(s,"gamma_data_Xud")
	Duu = ctn.load_as_np(s,"gamma_data_Duu")
	Ddd = ctn.load_as_np(s,"gamma_data_Ddd")
	Dud = ctn.load_as_np(s,"gamma_data_Dud")
	gamma_data_long_str_list[paid_order_ind,0,:,:] = Puu;
	gamma_data_long_str_list[paid_order_ind,1,:,:] = Pdd;
	gamma_data_long_str_list[paid_order_ind,2,:,:] = Pud;
	gamma_data_long_str_list[paid_order_ind,3,:,:] = Xud;
	gamma_data_long_str_list[paid_order_ind,4,:,:] = Duu;
	gamma_data_long_str_list[paid_order_ind,5,:,:] = Ddd;
	gamma_data_long_str_list[paid_order_ind,6,:,:] = Dud;
	comp_time_list[paid_order_ind] = ctn.load_as_np(s,"computation_time")[0,0]; 
	
	


#Determine data accuracy
diff_short_str=0.0
diff_long_str=0.0
for paid_order_ind in range(len(paid_order_list)):
	for paid_order_ind_2 in range(len(paid_order_list)):
		diff = np.amax(abs(gamma_data_short_str_list[paid_order_ind_2,:]-gamma_data_short_str_list[paid_order_ind,:]))
		print(diff)
		diff_short_str = max(diff_short_str,np.amax(abs(gamma_data_short_str_list[paid_order_ind_2,:]-gamma_data_short_str_list[paid_order_ind,:])))

		diff_long = np.amax(abs(gamma_data_long_str_list[paid_order_ind_2,:]-gamma_data_long_str_list[paid_order_ind,:]))
		print(diff_long)
		diff_long_str = max(diff_long_str,np.amax(abs(gamma_data_long_str_list[paid_order_ind_2,:]-gamma_data_long_str_list[paid_order_ind,:])))
print(diff_short_str)
print(diff_long_str)



fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(paid_order_list, comp_time_list, marker='o')
axis.set_xlabel(r'PAID order')
axis.set_ylabel(r'Computation_time/seconds')
axis.text(8,4500,r'QPC, 31 sites, $8$ Nodes')
plt.savefig('qpc_omp_standard_paid_order_scaling.pdf', format = 'pdf')
plt.show()
