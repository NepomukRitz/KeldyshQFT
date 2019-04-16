from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

s1="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_L1_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.500000_U1_0.000000_Xi_5.000000_T_0.002500_Lambda_ini_99999.500001_Lambda_fin_0.000000.mat"

A = sio.loadmat(s1)

L=A['L'][0,0]
N=A['N'][0,0]
pos_NfbX_0 = A['pos_NfbX_0'][0,0]

#Check correspondence static and central dynamic part:

#diff = amax(A['dy_new_data_long_str'][0,4]['m'][L,L]['m'] - A['dy_new_data_short_str'][0,4]['m'][0,pos_NfbX_0]['m'] ) 
static_center = A['gamma_data_long_str'][0,6]['m'][L,L]['m']; 
dynamic_center = A['gamma_data_short_str'][0,6]['m'][0,pos_NfbX_0]['m'];

diff = amax(absolute(static_center - dynamic_center))

print(diff)

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (3.375,2))

#axis.imshow(A['gamma_data_long_str'][0,0]['m'][L+1,L+1]['m'])
axis.imshow(A['gamma_data_long_str'][0,0]['m'][L+1,L+1]['m'], interpolation='None')
plt.show();
