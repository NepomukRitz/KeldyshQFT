import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

s="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi/dsfRG_mpi_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"
s_2_nodes="/gpfs/work/hmu26/hmu261/DATA/Test_full_mpi/dsfRG_mpi_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_2.mat"

performance_track = ctn.load_as_np(s,"performance_track");
performance_track_2_nodes = ctn.load_as_np(s_2_nodes,"performance_track");


fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))

axis[0].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,1], color='blue',marker='o')
axis[0].semilogx(performance_track_2_nodes[0,:,0,0],performance_track_2_nodes[0,:,0,1], color='blue',marker='D')
axis[0].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,4], color='green',marker='o')
axis[0].semilogx(performance_track_2_nodes[0,:,0,0],performance_track_2_nodes[0,:,0,4], color='green',marker='D')
axis[0].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,7], color='red',marker='o')
axis[0].semilogx(performance_track_2_nodes[0,:,0,0],performance_track_2_nodes[0,:,0,7], color='red',marker='D')
#
axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,2], color='orange',marker='o')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,2], color='blue',marker='o')
axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,5], color='orange',marker='o')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,5], color='green',marker='o')
axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,8], color='orange',marker='o')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,8], color='red',marker='o')

axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,3], color='orange',marker='o',linestyle='--')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,3], color='blue',marker='o',linestyle='--')
axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,6], color='orange',marker='o',linestyle='--')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,6], color='green',marker='o',linestyle='--')
axis[1].semilogx(performance_track_2_nodes[0,:,0,0],performance_track[0,:,0,9], color='orange',marker='o',linestyle='--')
axis[1].semilogx(performance_track[0,:,0,0],performance_track[0,:,0,9], color='red',marker='o',linestyle='--')
plt.show()

