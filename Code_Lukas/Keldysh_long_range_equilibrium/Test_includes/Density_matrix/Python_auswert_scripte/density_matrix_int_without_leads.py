from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec

#N=15
N=8
Nges=2*N+1
Vg=0.0
mu=0.0
U0=0.8
#T_list=[0.0,0.39810000,1.0]
T_list=[0.1]
#T_list=[0.015800, 0.020000, 0.039800000]
#T_list=[0.398100,0.631000, 1.000000]
num_freq_pre=600000
seed=1

H0=empty(2*N+1)
ldos_integrated_int_T=empty([len(T_list),2*N+1,2*N+1])
density_int_T=empty([len(T_list),2*N+1,2*N+1])
eig_density_int_T=empty([len(T_list),2*N+1])
for ind_T in range(len(T_list)):
 	T = T_list[ind_T]
	foldername="/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/Without_leads_reg3/"
	filename="density_matrix_N%d_Vg%f_mu%f_T%f_U0%f_num_freq_pre%d_seed%d.mat" % (N,Vg,mu,T,U0,num_freq_pre,seed)
	s=foldername + filename
	A = sio.loadmat(s)
	ldos_integrated_int_T[ind_T,:,:]=A['ldos_integrated']
	density_int_T[ind_T,:,:]=A['density']
	eig_density_int_T[ind_T,:]=A['eig_density'][0,:]
	H0 = A['hamiltonian']

ldos_integrated_non_int_T=empty([len(T_list),2*N+1,2*N+1])
density_non_int_T=empty([len(T_list),2*N+1,2*N+1])
eig_density_non_int_T=empty([len(T_list),2*N+1])
for ind_T in range(len(T_list)):
 	T = T_list[ind_T]
	foldername="/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/Without_leads/"
	filename="density_matrix_non_int_N%d_Vg%f_mu%f_T%f_num_freq_pre%d_seed%d.mat" % (N,Vg,mu,T,num_freq_pre,seed)
	s=foldername + filename
	A = sio.loadmat(s)
	ldos_integrated_non_int_T[ind_T,:,:]=A['ldos_integrated']
	density_non_int_T[ind_T,:,:]=A['density']
	eig_density_non_int_T[ind_T,:]=A['eig_density'][0,:]
sites=arange(-N,N+1,1)

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (12,8))
#fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))

axis[0].plot(sites,real(diag(H0)),color='black')
for ind_T in range(len(T_list)):
	axis[0].plot(sites,real(diag(ldos_integrated_int_T[ind_T,:,:])),linestyle='--')
axis[0].set_prop_cycle(None)
for ind_T in range(len(T_list)):
	axis[0].plot(sites,real(diag(density_int_T[ind_T,:,:])),label=T_list[ind_T])
legend = axis[0].legend(loc='upper left', fontsize=5)
plt.gca().set_prop_cycle(None)

for ind_T in range(len(T_list)):
 	axis[1].plot(arange(0,2*N+1,1),eig_density_int_T[ind_T,:],marker='o', label=T_list[ind_T])
 	#axis.plot(arange(0,2*N+1,1),eig_density_int_T[ind_T,:],marker='o', label=T_list[ind_T])
legend = axis[1].legend(loc='upper left', fontsize=5)
#legend = axis.legend(loc='upper left', fontsize=5)
plt.gca().set_prop_cycle(None)
for ind_T in range(len(T_list)):
 	axis[1].plot(arange(0,2*N+1,1),eig_density_non_int_T[ind_T,:],marker='x',linestyle=':')
 	#axis.plot(arange(0,2*N+1,1),eig_density_non_int_T[ind_T,:],marker='x',linestyle=':')
plt.gca().set_prop_cycle(None)


savename='density_matrix_int_without_leads_Vg%f_seed%d.pdf' % (Vg,seed)
plt.savefig(savename, format = 'pdf')
plt.show();
