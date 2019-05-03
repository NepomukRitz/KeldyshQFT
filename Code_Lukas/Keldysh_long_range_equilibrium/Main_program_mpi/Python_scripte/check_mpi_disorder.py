import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

	
s_folder="/gpfs/work/hmu26/hmu261/DATA/Paid_systematic_test_correct/N15_Disorder/Standard_code/"
s1_file="dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_1.200000_h_0.000000_mu_0.000000_T_0.100000_U0_0.800000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_8.mat"

s=s_folder + s1_file

H0 = np.real(ctn.load_as_np(s,"hamiltonian"));

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(np.diag(H0))
plt.show()

