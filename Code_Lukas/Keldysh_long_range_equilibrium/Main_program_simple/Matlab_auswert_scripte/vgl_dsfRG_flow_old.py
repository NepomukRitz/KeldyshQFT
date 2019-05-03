from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (3.375,2))
		
s1 = "/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/flow_old_L0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.500000_T_0.030000_Lambda_ini_9999.500008_Lambda_fin_0.000000.mat"

s2 = "/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_L0_Lu0_N4_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.500000_U1_0.000000_Xi_5.000000_T_0.030000_Lambda_ini_9999.500008_Lambda_fin_0.000000.mat"


A1 = sio.loadmat(s1)
A2 = sio.loadmat(s2)

wbP_1 = A1['wbP'][0,:]
wbP_2 = A2['wbP'][0,:]
wbX_1 = A1['wbX'][0,:]
wbX_2 = A2['wbX'][0,:]
wf_1 = A1['wf'][0,:]
wf_2 = A2['wf'][0,:]
diff_wbP = amax(absolute( wbP_1 - wbP_2 ))
diff_wbX = amax(absolute( wbX_1 - wbX_2 ))
Nff = A1['Nff'][0,0]
N = A1['N'][0,0]
NfbP_1 = A1['NfbP'][0,0]
NfbP_2 = A2['NfbP'][0,0]
NfbX_1 = A1['NfbX'][0,0]
NfbX_2 = A2['NfbX'][0,0]
pos_NfbP_2mu_1 = A1['pos_NfbP_2mu'][0,0]
pos_NfbP_2mu_2 = A2['pos_NfbP_2mu'][0,0]
pos_NfbX_0_1 = A1['pos_NfbX_0'][0,0]
pos_NfbP_0_2 = A2['pos_NfbX_0'][0,0]
Nges = A1['Nges'][0,0]
print(diff_wbP)
print(diff_wbX)

a_p_1 = empty((NfbP_1,Nges,Nges),dtype=complex64)
a_p_2 = empty((NfbP_2,Nges,Nges),dtype=complex64)
a_x_1 = empty((NfbX_1,Nges,Nges),dtype=complex64)
a_x_2 = empty((NfbX_2,Nges,Nges),dtype=complex64)
a_d_1 = empty((NfbX_1,Nges,Nges),dtype=complex64)
a_d_2 = empty((NfbX_2,Nges,Nges),dtype=complex64)
self_1 = empty((Nff,Nges,Nges),dtype=complex64)
self_2 = empty((Nff,Nges,Nges),dtype=complex64)

for indf in range(NfbP_1):
	a_p_1[indf,:,:] = A1['y_old'][0,2]['m'][0,indf]['m'][:,:]-0.25*A1['U'][:,:]
for indf in range(NfbP_2):
	a_p_2[indf,:,:] = A2['gamma_data_short_str'][0,2]['m'][0,indf]['m'][:,:]

for indf in range(NfbX_1):
	a_x_1[indf,:,:] = A1['y_old'][0,3]['m'][0,indf]['m'][:,:]-0.25*A1['U'][:,:]
for indf in range(NfbX_2):
	a_x_2[indf,:,:] = A2['gamma_data_short_str'][0,3]['m'][0,indf]['m'][:,:]

for indf in range(NfbX_1):
	a_d_1[indf,:,:] = A1['y_old'][0,4]['m'][0,indf]['m'][:,:]
for indf in range(NfbX_2):
	a_d_2[indf,:,:] = A2['gamma_data_short_str'][0,4]['m'][0,indf]['m'][:,:]

for indf in range(Nff):
	self_1[indf,:,:] = A1['y_old'][0,0]['m'][0,indf]['m'][:,:]
for indf in range(NfbX_2):
	self_2[indf,:,:] = A2['gamma_data_short_str'][0,7]['m'][0,indf]['m'][:,:]

diff_a_p = amax(absolute( a_p_1 - a_p_2))
diff_a_x = amax(absolute( a_x_1 - a_x_2))
diff_a_d = amax(absolute( a_d_1 - a_d_2))
diff_self = amax(absolute( self_1 - self_2))
print("diff_a_p={0:.9f}".format(diff_a_p))
print("diff_a_x={0:.9f}".format(diff_a_x))
print("diff_a_d={0:.9f}".format(diff_a_d))
print("diff_self={0:.9f}".format(diff_self))



#axis.plot(wbP_1,imag(a_p_1[:,N,N]), color='blue')
#axis.plot(wbP_2,imag(a_p_2[:,N,N]), color='red')
#axis.plot(wbP_1,real(a_p_1[:,N,N]), color='blue')
#axis.plot(wbP_2,real(a_p_2[:,N,N]), color='red')

#axis.plot(wbX_1,imag(a_x_1[:,N,N]), color='blue')
#axis.plot(wbX_2,imag(a_x_2[:,N,N]), color='red')
#axis.plot(wbX_1,real(a_x_1[:,N,N]), color='blue')
#axis.plot(wbX_2,real(a_x_2[:,N,N]), color='red')

#axis.plot(wbX_1,imag(a_d_1[:,N,N]), color='blue')
#axis.plot(wbX_2,imag(a_d_2[:,N,N]), color='red')
#axis.plot(wbX_1,real(a_d_1[:,N,N]), color='blue')
#axis.plot(wbX_2,real(a_d_2[:,N,N]), color='red')

axis.plot(wf_1,imag(self_1[:,N,N]), color='blue')
axis.plot(wf_2,imag(self_2[:,N,N]), color='red')
axis.plot(wf_1,real(self_1[:,N,N]), color='blue')
axis.plot(wf_2,real(self_2[:,N,N]), color='red')

axis.set_xlim([-5,5])
#
#plt.savefig('vgl_dsfRG_flow_old.pdf', format = 'pdf')
plt.show();
