import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_freq_str as Ex_freq_str 

s="/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Test_python_files/ex_compute_stat_bubble.mat"
[L, N, L_structure, wb, bubble_dyn, bubble_stat] = Ex_freq_str.load(s,"P_bubble")
print("np.shape(bubble_dyn)=", np.shape(bubble_dyn))
print("np.shape(bubble_stat)=", np.shape(bubble_stat))

tmp = ctn.load_as_np("ex_compute_stat_bubble.mat","pos_NfbP_2mu")
pos_NfbP_2mu = int(np.real(tmp[0,0]))

#print("bubble_stat=", bubble_stat)
#print("bubble_dyn[0]=",bubble_dyn[0])
Li = L_structure[pos_NfbP_2mu]
print("Li=", Li)
print(np.amax(abs(bubble_dyn[pos_NfbP_2mu][Li,Li] - bubble_stat[L,L])))

