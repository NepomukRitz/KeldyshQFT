import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str
import One_particle_conductance_to_numpy as Opctn 


s="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/X_L0_Lu0_N30_Nff1500_NfbP125_NfbX125_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4500_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"

H = np.real(ctn.load_as_np(s,"hamiltonian"))
mu_list = np.linspace(-1.6,-1.4,20)
h=0.0;
cond_up_list = Opctn.get_cond_curve_non_int(1,mu_list,h,H)
cond_down_list = Opctn.get_cond_curve_non_int(0,mu_list,h,H)

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))

axis.plot(mu_list,cond_up_list+cond_down_list,color='blue',marker='o')

#axis.set_xlim([-5,5])

plt.show()


