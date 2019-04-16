import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

number_of_nodes = 1 

s2_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/Standard_code/"
s2_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (number_of_nodes)
s2=s2_folder + s2_file
performance_track_one_node = ctn.load_as_np(s2,"performance_track");

number_of_nodes = 8 
s2_folder="/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_QPC/Standard_code/"
s2_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_%d.mat" % (number_of_nodes)
s2=s2_folder + s2_file
performance_track_eight_node = ctn.load_as_np(s2,"performance_track");


fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))

axis[0].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,1], color='blue',marker='o')
axis[0].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,4], color='green',marker='o')
axis[0].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,7], color='red',marker='o')

axis[0].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,1], color='blue',marker='D')
axis[0].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,4], color='green',marker='D')
axis[0].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,7], color='red',marker='D')

axis[1].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,0] -performance_track_eight_node[0,:,0,0] , color='blue',marker='o')
#axis[1].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,5]*31*31, color='green',marker='o')
#axis[1].semilogx(performance_track_one_node[0,:,0,0],performance_track_one_node[0,:,0,8]*31*31, color='red',marker='o')

#axis[1].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,2]*31*31, color='lightblue',marker='D',linestyle='--')
#axis[1].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,5]*31*31, color='green',marker='D')
#axis[1].semilogx(performance_track_eight_node[0,:,0,0],performance_track_eight_node[0,:,0,8]*31*31, color='red',marker='D')

plt.show()

