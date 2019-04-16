from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

s1="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_evaluate_right_hand_side.mat"

A = sio.loadmat(s1)
L=A['L'][0,0]
N=A['N'][0,0]
pos_NfbX_0 = A['pos_NfbX_0'][0,0]

#Check correspondence static and central dynamic part:

#diff = amax(A['dy_new_data_long_str'][0,4]['m'][L,L]['m'] - A['dy_new_data_short_str'][0,4]['m'][0,pos_NfbX_0]['m'] ) 
static_center = A['dy_new_data_long_str'][0,6]['m'][L,L]['m']; 
dynamic_center = A['dy_new_data_short_str'][0,6]['m'][0,pos_NfbX_0]['m'];

diff = amax(absolute(static_center - dynamic_center))

print(diff)
