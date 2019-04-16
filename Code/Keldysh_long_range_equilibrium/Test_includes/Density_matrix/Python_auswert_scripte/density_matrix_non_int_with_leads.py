from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec

N=15
Vg=1.2
mu=0.0
T_list=[0.0,0.5,1.0]
num_freq_pre=600000

H0=empty(2*N+1)
ldos_integrated_T=empty([len(T_list),2*N+1,2*N+1])
density_T=empty([len(T_list),2*N+1,2*N+1])
eig_density_T=empty([len(T_list),2*N+1])
for ind_T in range(len(T_list)):
 	T = T_list[ind_T]
	foldername="/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/With_leads/"
	filename="density_matrix_non_int_N%d_Vg%f_mu%f_T%f_num_freq_pre%d.mat" % (N,Vg,mu,T,num_freq_pre)
	s=foldername + filename
	A = sio.loadmat(s)
	ldos_integrated_T[ind_T,:,:]=A['ldos_integrated']
	density_T[ind_T,:,:]=A['density']
	eig_density_T[ind_T,:]=A['eig_density'][0,:]
	H0 = A['H0']
sites=arange(-N,N+1,1)

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (3.375,2))

axis[0].plot(sites,real(diag(H0)),color='black')
axis[0].plot(sites,real(diag(ldos_integrated_T[0,:,:])),color='blue',linestyle='--')
axis[0].plot(sites,real(diag(density_T[0,:,:])),color='blue')
axis[0].plot(sites,real(diag(ldos_integrated_T[1,:,:])),color='green',linestyle='--')
axis[0].plot(sites,real(diag(density_T[1,:,:])),color='green')
axis[0].plot(sites,real(diag(ldos_integrated_T[2,:,:])),color='red',linestyle='--')
axis[0].plot(sites,real(diag(density_T[2,:,:])),color='red')

for ind_T in range(len(T_list)):
 	axis[1].plot(arange(0,2*N+1,1),eig_density_T[ind_T,:],marker='o')
plt.savefig('density_matrix_non_int_with_leads.pdf', format = 'pdf')
plt.show();
