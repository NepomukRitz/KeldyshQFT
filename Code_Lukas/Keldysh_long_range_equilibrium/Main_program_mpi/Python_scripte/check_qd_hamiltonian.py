import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt

L_list=[0]
#nodes=[1,5,10,20,30,40,50,60,70]
nodes=[5]

s="/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/Schar_dsfRG_mpi_L0_Lu0_N40_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg2.0_h0.0_T0.0025_U00.0_U10.0_Xi5.0_Vsg_pot_width10_mu0.0_0.01_0.0_nodes10/dsfRG_mpi_L0_Lu0_N40_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_2.000000_h_0.000000_mu_0.000000_T_0.002500_U0_0.000000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_10.mat"

s2="/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/Schar_dsfRG_mpi_L0_Lu0_N40_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg-1.0_h0.0_T0.0025_U00.0_U10.0_Xi5.0_Vsg_pot_width10_mu-1.0_0.01_-1.0_nodes10/dsfRG_mpi_L0_Lu0_N40_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_-1.000000_h_0.000000_mu_-1.000000_T_0.002500_U0_0.000000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_10.mat"

s3="/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/Schar_dsfRG_mpi_L0_Lu0_N40_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.9_h0.0_T0.0025_U00.0_U10.0_Xi5.0_Vsg_pot_width10_mu-1.0_0.01_-1.0_nodes10/dsfRG_mpi_L0_Lu0_N40_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.900000_h_0.000000_mu_-1.000000_T_0.002500_U0_0.000000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_10.mat"

s4="/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/Schar_dsfRG_mpi_L0_Lu0_N40_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.9_h0.0_T0.0025_U00.0_U10.0_Xi5.0_Vsg_pot_width20_mu-1.0_0.01_-1.0_nodes10/dsfRG_mpi_L0_Lu0_N40_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.900000_h_0.000000_mu_-1.000000_T_0.002500_U0_0.000000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_10.mat"

s5="/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/Schar_dsfRG_mpi_L0_Lu0_N40_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_mu-1.0_h0.0_T0.0025_U00.0_U10.0_Xi5.0_Vsg_pot_width20_Vg0.6_0.01_0.6_nodes5/dsfRG_mpi_L0_Lu0_N40_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.600000_h_0.000000_mu_-1.000000_T_0.002500_U0_0.000000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_5.mat"

H0 = ctn.load_as_np(s5,"hamiltonian")

print(H0)

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(np.real(np.diag(H0,1)))

plt.show()
