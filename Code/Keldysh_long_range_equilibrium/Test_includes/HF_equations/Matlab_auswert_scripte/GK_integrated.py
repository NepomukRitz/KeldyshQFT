from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (3.375,2))

s1="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/HF_equations/hf_equations_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.500000_U1_0.000000_Xi_5.000000_T_0.030000_Lambda_0.000000010_acc_hf0.000100000.mat"
A1 = sio.loadmat(s1)

print(amax(absolute(A1['Selfenergy'])))
#axis[0].imshow(real(A1['Selfenergy']),interpolation='None')

#plt.savefig('GK_integrated.pdf', format = 'pdf')
#plt.show()
