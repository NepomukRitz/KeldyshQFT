import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import os

def read(sdir):
	#directory_in_str="/gpfs/work/hmu26/hmu261/DATA/Unit_tests/Conductance_lead_system/Auto/"
	directory = os.fsencode(sdir)
	
	mu_list=[];
	cond_list=[];
	z=0
	for file in os.listdir(directory):
		s = sdir + os.fsdecode(file)
		if s.endswith(".mat"): 
			mu = ctn.load_as_np(s,"mu")[0,0]
			mu_list.append(mu)
			cond = ctn.load_as_np(s,"cond_prov")[0,0]
			cond_list.append(cond)
			z=z+1;
	
	print(len(mu_list))
	
	mu_vector = np.empty([len(mu_list)])
	mu_vector, cond_vector = zip(*sorted(zip(mu_list, cond_list)))
	print(mu_vector)
	print(cond_vector)
	return [mu_vector,cond_vector]








fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
#
[mu_vector,cond_vector] = read("/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0_conductance/Conductance_schar_dsfRG_mpi_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/");
axis.plot(np.real(mu_vector),np.real(cond_vector), color='blue', marker='o')
[mu_vector,cond_vector] = read("/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0_conductance/Conductance_schar_dsfRG_mpi_L1_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/");
axis.plot(np.real(mu_vector),np.real(cond_vector), color='green', marker='o')
[mu_vector,cond_vector] = read("/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0_conductance/Conductance_schar_dsfRG_mpi_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/");
axis.plot(np.real(mu_vector),np.real(cond_vector), color='red', marker='o')
[mu_vector,cond_vector] = read("/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0_conductance/Conductance_schar_dsfRG_mpi_L10_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/");
axis.plot(np.real(mu_vector),np.real(cond_vector), color='orange', marker='o')
#[mu_vector,cond_vector] = read("/gpfs/work/hmu26/hmu261/DATA/Unit_tests/Conductance_lead_system/Auto/");
#axis[0].plot(np.real(mu_vector),np.real(cond_vector), color='blue', marker='o')
#axis[0].plot(freq[0,:],np.real(Integrated_vertex_contribution[0,:,15,15]), color='red', marker='o')
#
plt.savefig("qpc_various_L.pdf", format = 'pdf')
plt.show()

