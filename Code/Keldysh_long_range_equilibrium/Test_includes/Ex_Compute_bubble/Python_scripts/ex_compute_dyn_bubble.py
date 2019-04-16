import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import ex_str as ex

s="/p/home/jusers/weidinger1/juwels/fRG_development/Keldysh_long_range_equilibrium/Test_includes/Ex_Compute_bubble/Start_scripts/ex_compute_dyn_bubble.mat"
#[L, N, L_structure, wb, bubble_dyn, bubble_stat] = ex.load_ex_freq_str(s,"P_bubble")
#[L, N, L_structure, wb, bubble_dyn, bubble_stat] = ex.load_ex_freq_str(s,"X_bubble")
[L, N, L_structure, wb, bubble_dyn, bubble_stat] = ex.load_ex_freq_str(s,"Ddd_bubble")
#print(L)
#print(N)
#print(L_structure)
#print(bubble_dyn)
#print(bubble_stat)
#print(bubble_dyn[:,0,0,N,N])

#bubble_old = ctn.load_as_np(s,"P_Bubble_bubEq_precomputed")
#bubble_old = ctn.load_as_np(s,"X_Bubble_bubEq_precomputed")
bubble_old = ctn.load_as_np(s,"Ddd_Bubble_bubEq_precomputed")
#print(np.shape(bubble_dyn))
#print(np.shape(wb))
fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (15,10))
axis.plot(wb,np.imag(bubble_dyn[:,0,0,N,N]), color='blue', marker='o')
axis.plot(wb,np.real(bubble_dyn[:,0,0,N,N]), color='blue', marker='o')
axis.plot(wb,np.imag(bubble_old[0,:,N,N]), color='red', marker='o')
axis.plot(wb,np.real(bubble_old[0,:,N,N]), color='red', marker='o')
axis.set_xlim([-5,5])
plt.show()
