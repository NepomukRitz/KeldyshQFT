import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

#s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test/N15_QPC/OMP_standard/"
s_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/OMP_standard/"

paid_order = 4
node_list=[1,2,4,8,16,32,64]

#Preload to determine dimensions:
s_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_%d.mat" % (paid_order,node_list[0])
s=s_folder + s_file
Nges = ctn.load_as_np(s,"Nges")[0,0].astype(int);
L = ctn.load_as_np(s,"L")[0,0];
wf = ctn.load_as_np(s,"wf")[0,:];
Puu = ctn.load_as_np(s,"gamma_data_Puu")
Block_dim = np.shape(Puu)[0]
Nff_real = len(wf)

gamma_data_short_str_list=np.empty([len(node_list),1,9, 1, Nff_real,Nges,Nges],dtype=complex)
gamma_data_long_str_list=np.empty([len(node_list),7,Block_dim,Block_dim])

comp_time_list =np.empty(len(node_list))

#Load all data:
for node_ind in range(len(node_list)):
	node_number = node_list[node_ind]
	s_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_%d.mat" % (paid_order,node_number)
	s=s_folder + s_file
	gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
	gamma_data_short_str_list[node_ind,:,:,:,:,:,:] = gamma_data_short_str[:,:,:,:,:,:]
	Puu = ctn.load_as_np(s,"gamma_data_Puu")
	Pdd = ctn.load_as_np(s,"gamma_data_Pdd")
	Pud = ctn.load_as_np(s,"gamma_data_Pud")
	Xud = ctn.load_as_np(s,"gamma_data_Xud")
	Duu = ctn.load_as_np(s,"gamma_data_Duu")
	Ddd = ctn.load_as_np(s,"gamma_data_Ddd")
	Dud = ctn.load_as_np(s,"gamma_data_Dud")
	gamma_data_long_str_list[node_ind,0,:,:] = Puu;
	gamma_data_long_str_list[node_ind,1,:,:] = Pdd;
	gamma_data_long_str_list[node_ind,2,:,:] = Pud;
	gamma_data_long_str_list[node_ind,3,:,:] = Xud;
	gamma_data_long_str_list[node_ind,4,:,:] = Duu;
	gamma_data_long_str_list[node_ind,5,:,:] = Ddd;
	gamma_data_long_str_list[node_ind,6,:,:] = Dud;
	comp_time_list[node_ind] = ctn.load_as_np(s,"computation_time")[0,0]; 

#Determine data accuracy
diff_short_str=0.0
diff_long_str=0.0
for node_ind in range(len(node_list)):
	for node_ind_2 in range(len(node_list)):
		diff = np.amax(abs(gamma_data_short_str_list[node_ind_2,:]-gamma_data_short_str_list[node_ind,:]))
		print(diff)
		diff_short_str = max(diff_short_str,np.amax(abs(gamma_data_short_str_list[node_ind_2,:]-gamma_data_short_str_list[node_ind,:])))

		diff_long = np.amax(abs(gamma_data_long_str_list[node_ind_2,:]-gamma_data_long_str_list[node_ind,:]))
		print(diff_long)
		diff_long_str = max(diff_long_str,np.amax(abs(gamma_data_long_str_list[node_ind_2,:]-gamma_data_long_str_list[node_ind,:])))
print(diff_short_str)
print(diff_long_str)

#Load standard code data
gamma_data_short_str_standard_list=np.empty([len(node_list),1,9, 1, Nff_real,Nges,Nges],dtype=complex)
gamma_data_long_str_standard_list=np.empty([len(node_list),7,Block_dim,Block_dim])

comp_time_standard_list =np.empty(len(node_list))
for node_ind in range(len(node_list)):
	node_number = node_list[node_ind]
	#s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test/N15_QPC/Standard_code/"
	s_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/Standard_code/"
	s_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (node_number)
	s=s_folder + s_file
	gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
	gamma_data_short_str_standard_list[node_ind,:,:,:,:,:,:] = gamma_data_short_str[:,:,:,:,:,:]
	Puu = ctn.load_as_np(s,"gamma_data_Puu")
	Pdd = ctn.load_as_np(s,"gamma_data_Pdd")
	Pud = ctn.load_as_np(s,"gamma_data_Pud")
	Xud = ctn.load_as_np(s,"gamma_data_Xud")
	Duu = ctn.load_as_np(s,"gamma_data_Duu")
	Ddd = ctn.load_as_np(s,"gamma_data_Ddd")
	Dud = ctn.load_as_np(s,"gamma_data_Dud")
	gamma_data_long_str_standard_list[node_ind,0,:,:] = Puu;
	gamma_data_long_str_standard_list[node_ind,1,:,:] = Pdd;
	gamma_data_long_str_standard_list[node_ind,2,:,:] = Pud;
	gamma_data_long_str_standard_list[node_ind,3,:,:] = Xud;
	gamma_data_long_str_standard_list[node_ind,4,:,:] = Duu;
	gamma_data_long_str_standard_list[node_ind,5,:,:] = Ddd;
	gamma_data_long_str_standard_list[node_ind,6,:,:] = Dud;
	comp_time_standard_list[node_ind] = ctn.load_as_np(s,"computation_time")[0,0]; 

#Determine data accuracy within standard_code
print("standard code:")
print(np.amax(abs(gamma_data_short_str_standard_list[0,0,3,:])))
diff_short_str=0.0
diff_long_str=0.0
for node_ind in range(len(node_list)):
	for node_ind_2 in range(len(node_list)):
		diff = np.amax(abs(gamma_data_short_str_standard_list[node_ind_2,:]-gamma_data_short_str_standard_list[node_ind,:]))
		print(node_ind)
		print(node_ind_2)
		print(diff)
		diff_short_str = max(diff_short_str,np.amax(abs(gamma_data_short_str_standard_list[node_ind_2,:]-gamma_data_short_str_standard_list[node_ind,:])))

		diff_long = np.amax(abs(gamma_data_long_str_standard_list[node_ind_2,:]-gamma_data_long_str_standard_list[node_ind,:]))
		print(diff_long)
		diff_long_str = max(diff_long_str,np.amax(abs(gamma_data_long_str_standard_list[node_ind_2,:]-gamma_data_long_str_standard_list[node_ind,:])))
print(diff_short_str)
print(diff_long_str)

#Compare data accuracy with PAID

print("comparison:")
diff = np.amax(abs(gamma_data_short_str_list[0,:]-gamma_data_short_str_standard_list[0,:]))
print(diff)

#Load OMP_inside data:

node_omp_inside_list=[4,8,16]

gamma_data_short_str_omp_inside_list=np.empty([len(node_omp_inside_list),1,9, 1, Nff_real,Nges,Nges],dtype=complex)
gamma_data_long_str_omp_inside_list=np.empty([len(node_omp_inside_list),7,Block_dim,Block_dim])

comp_time_omp_inside_list =np.empty(len(node_omp_inside_list))

for node_ind in range(len(node_omp_inside_list)):
	node_number = node_omp_inside_list[node_ind]
	s_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/OMP_inside/"
	s_file="dsfRG_mpi_with_paid_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_acc_p1.00e-04_acc_x1.00e-04_acc_s1.00e-04_acc_f1.00e-06_pd_order%d_number_of_nodes_%d.mat" % (paid_order,node_number)
	s=s_folder + s_file
	gamma_data_short_str = ctn.load_as_np(s,"gamma_data_short_str")
	gamma_data_short_str_omp_inside_list[node_ind,:,:,:,:,:,:] = gamma_data_short_str[:,:,:,:,:,:]
	Puu = ctn.load_as_np(s,"gamma_data_Puu")
	Pdd = ctn.load_as_np(s,"gamma_data_Pdd")
	Pud = ctn.load_as_np(s,"gamma_data_Pud")
	Xud = ctn.load_as_np(s,"gamma_data_Xud")
	Duu = ctn.load_as_np(s,"gamma_data_Duu")
	Ddd = ctn.load_as_np(s,"gamma_data_Ddd")
	Dud = ctn.load_as_np(s,"gamma_data_Dud")
	gamma_data_long_str_omp_inside_list[node_ind,0,:,:] = Puu;
	gamma_data_long_str_omp_inside_list[node_ind,1,:,:] = Pdd;
	gamma_data_long_str_omp_inside_list[node_ind,2,:,:] = Pud;
	gamma_data_long_str_omp_inside_list[node_ind,3,:,:] = Xud;
	gamma_data_long_str_omp_inside_list[node_ind,4,:,:] = Duu;
	gamma_data_long_str_omp_inside_list[node_ind,5,:,:] = Ddd;
	gamma_data_long_str_omp_inside_list[node_ind,6,:,:] = Dud;
	comp_time_omp_inside_list[node_ind] = ctn.load_as_np(s,"computation_time")[0,0]; 

#Compare data accuracy between inside and normal OMP mode 

print("comparison inside:")
diff = np.amax(abs(gamma_data_short_str_list[0,:]-gamma_data_short_str_omp_inside_list[0,:]))
print(diff)

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(node_list, comp_time_standard_list, marker='o', color='blue')
axis.plot(node_list, comp_time_list, marker='o',color='green')
axis.plot(node_omp_inside_list, comp_time_omp_inside_list, marker='o',color='red')
axis.set_xlabel(r'Node number')
axis.set_ylabel(r'Computation_time/seconds')
axis.text(30,12000,r'QPC, $31$ sites')
axis.text(15,6000,r'PAID: OMP inside')
axis.text(20,1600,r'PAID: OMP outside')
axis.text(10,10,r'Old integrator')
plt.savefig('qpc_omp_standard_node_scaling.pdf', format = 'pdf')
plt.show()
