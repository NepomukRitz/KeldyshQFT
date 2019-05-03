import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import ex_str as ex

s="/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Test_includes/Ex_Compute_bubble/Start_scripts/ex_compute_stat_bubble.mat"
[L, N, L_structure, wb, bubble_dyn, bubble_stat] = ex.load_ex_freq_str(s,"P_bubble")

tmp = ctn.load_as_np(s,"pos_NfbP_2mu")
pos_NfbP_2mu = int(tmp[0,0])

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.imshow(np.real(bubble_stat[L,L]))
plt.show()
