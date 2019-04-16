from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec

#s="/naslx/projects/uh3o1/ri26yad/DATA/Keldysh_short_range_non_eq_disorder/ChainNEMoreFreqclean_disorder_central_build_N31_Vg1.200000_muL0.00000_muR0.00000_TL0.01500_TR0.01500_h0.00000_U0.800000_seed5.mat"	
#s="/naslx/projects/uh3o1/ri26yad/DATA/Keldysh_short_range_non_eq_disorder/ChainNEMoreFreqclean_disorder_central_build_N31_Vg0.250000_muL-1.47500_muR-1.47500_TL0.10000_TR0.10000_h0.00000_U0.800000_seed9999.mat"	
s="/naslx/projects/uh3o1/ri26yad/DATA/NE/Vg_25/ChainNEMoreFreqclean_N31_Vg0.250000_muL-1.47500_muR-1.47500_TL0.01000_TR0.01000_h0.00000_U0.550000.mat"

s_eq="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Conductance_zero_mag/conductance_zero_mag_N15_Vg0.250000_mu-1.475000_T0.010000_U00.550000.mat"	


A = sio.loadmat(s)

Nges = A['N'][0,0]
N=(Nges-1)/2
bP = A['bP']
wbP = A['wbP'][0,:]

p=empty([len(wbP)], dtype=complex)
for ind in range(len(wbP)):
 	p[ind] = bP[0,ind]['m'][N,N]

B = sio.loadmat(s_eq)
y= B['y']
wf_eq = B['wf'][0,:]
wbP_eq = B['wbP'][0,:]
wbX_eq = B['wbX'][0,:]
N_eq = B['N'][0,0]

p_eq=empty([len(wbP_eq)], dtype=complex)
for ind in range(len(wbP_eq)):
 	p_eq[ind] = y[0,5]['m'][0,ind]['m'][N_eq,N_eq]



plt.plot(wbP,imag(p),color='red',marker='o')
plt.plot(wbP_eq,imag(p_eq),color='blue',marker='o')
plt.xlim([-5, 5])
plt.show();
