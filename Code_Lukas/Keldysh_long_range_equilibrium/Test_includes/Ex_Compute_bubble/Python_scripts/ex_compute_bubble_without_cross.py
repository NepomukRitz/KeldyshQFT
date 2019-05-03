import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import ex_str as ex

s="/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Test_includes/Ex_Compute_bubble/Start_scripts/ex_compute_bubble_without_cross.mat"
[L, N, L_structure, wb, bubble_dyn, bubble_stat] = ex.load_ex_freq_str(s,"P_bubble")
[L_wc, N_wc, L_structure_wc, wb, bubble_dyn_wc, bubble_stat_wc] = ex.load_ex_freq_str(s,"P_bubble_wc")

tmp = ctn.load_as_np(s,"pos_NfbP_2mu")
pos_NfbP_2mu = int(np.real(tmp[0,0]))
pos_error = 50


fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))
axis[0].imshow(np.real(bubble_dyn[pos_error][1,3]))
axis[1].imshow(np.real(bubble_dyn_wc[pos_error][1,3]))
plt.show()
