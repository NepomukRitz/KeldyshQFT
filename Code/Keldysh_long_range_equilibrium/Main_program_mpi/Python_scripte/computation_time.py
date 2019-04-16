import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

L_list=[0,1,5]
#nodes=[1,5,10,20,30,40,50,60,70]
nodes=[5,10,20]


comp_time = np.empty([len(L_list),len(nodes)])
for ind_l in range(len(L_list)):
	for ind_n in range(len(nodes)):
		DIR_DATA="/work/hmu26/hmu261/DATA/Static_long_range/QPC/Test_parallelization_nodes/%d_nodes/Schar_dsfRG_mpi_L%d_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.03_U00.8_U10.0_Xi5.0_mu-1.475_0.01_-1.475_nodes%d/" % (nodes[ind_n],L_list[ind_l],nodes[ind_n])
		FILENAME="dsfRG_mpi_L%d_Lu0_N30_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.030000_U0_0.800000_U1_0.000000_Xi_5.000000_Lambda_ini_9999.500008_Lambda_fin_0.000000_number_of_nodes_%d.mat" %(L_list[ind_l],nodes[ind_n])
		s=DIR_DATA+FILENAME
		print(s)
		computation_time = ctn.load_as_np(s,"computation_time")
		comp_time[ind_l,ind_n] = np.real(ctn.load_as_np(s,"computation_time"))


fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(nodes,comp_time[0,:], marker='o',color='blue')
axis.plot(nodes,comp_time[1,:], marker='o',color='red')
axis.plot(nodes,comp_time[2,:], marker='o',color='green')
plt.show()
