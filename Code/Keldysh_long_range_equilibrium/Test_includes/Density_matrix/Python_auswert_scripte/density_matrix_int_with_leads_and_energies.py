from numpy import *
import scipy.io as sio
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.gridspec as gridspec

original_data_foldername="/naslx/projects/uh3o1/ri26yad/DATA/Keldysh_long_range_equilibrium/Full_select_data/"

def fermi(f,mu,T):
	return 1.0/(1.0 + exp((f-mu)/T))

#N=15
N=30
Vg=1.2
mu=0.0
U0=0.8
T_list=[1.0]
T_list=[0.015800, 0.020000, 0.025100, 0.039800000,0.100000, 0.158500, 0.251200,0.398100,0.631000, 1.000000]
#T_list=[0.015800, 0.398100,0.631000, 1.000000]
#T_list=[0.015800]
num_freq_pre=600000

H0=empty(2*N+1)
ldos_integrated_int_T=empty([len(T_list),2*N+1,2*N+1])
density_int_T=empty([len(T_list),2*N+1,2*N+1])
eig_density_int_T=empty([len(T_list),2*N+1])
associated_energies_int_T=empty([len(T_list),2*N+1])
associated_energies_int_T_tilde=empty([len(T_list),2*N+1])
for ind_T in range(len(T_list)):
 	print(ind_T)
 	T = T_list[ind_T]
	foldername="/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/With_leads/"
	filename="density_matrix_N%d_Vg%f_mu%f_T%f_U0%f_num_freq_pre%d.mat" % (N,Vg,mu,T,U0,num_freq_pre)
	s=foldername + filename
	A = sio.loadmat(s)
	ldos_integrated_int_T[ind_T,:,:]=A['ldos_integrated']
	density_int_T[ind_T,:,:]=A['density']
	eig_density_int_T[ind_T,:]=A['eig_density'][0,:]
	H0 = A['hamiltonian']
	w, v = linalg.eigh(density_int_T[ind_T,:,:])
	associated_energies = real(diag(matmul(matmul(v.transpose(),H0),v)))
	associated_energies_int_T[ind_T,:] = associated_energies 
	
	original_data_filename="Chain_disorder_central_build_N61_Vg1.200000_Vsg0.000000_width0_mu0.000000000_T%.9f_h0.000000000_U 0.800000_seed5.mat" % (T_list[ind_T])
	s2 = original_data_foldername + original_data_filename

	B = sio.loadmat(s2)
	wf = real(B['wf'][0,:])
	Sigma = empty([len(wf),2*N+1,2*N+1],dtype=complex)
	Sigma_associated = empty([len(associated_energies),2*N+1,2*N+1], dtype=complex)
	for ind_f in range(len(wf)):
	 	for ind_i in range(2*N+1):
	 		for ind_j in range(2*N+1):
				Sigma[ind_f,ind_i,ind_j] = B['m'][0,0]['m'][0,ind_f]['m'][ind_i,ind_j] 		

	for ind_as in range(len(associated_energies)):
		for ind_f in range(len(wf)):
			if (wf[ind_f] <= associated_energies[ind_as]) and (wf[ind_f+1] >= associated_energies[ind_as]) :
			 	Sigma_associated[ind_as,:,:] = Sigma[ind_f,:,:] 
		
	for ind_as in range(len(associated_energies)):
		associated_energies_int_T_tilde[ind_T,ind_as] = matmul(matmul(v[:,ind_as].transpose(),Sigma_associated[ind_as,:,:]),v[:,ind_as]) 
		

	
wf_fermi=arange(-5.,5.,0.01)		
fermi_val = empty([len(T_list),len(wf_fermi)])
	
ldos_integrated_non_int_T=empty([len(T_list),2*N+1,2*N+1])
density_non_int_T=empty([len(T_list),2*N+1,2*N+1])
eig_density_non_int_T=empty([len(T_list),2*N+1])
associated_energies_non_int_T=empty([len(T_list),2*N+1])
for ind_T in range(len(T_list)):
 	T = T_list[ind_T]
	print(ind_T)
	foldername="/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/With_leads/"
	filename="density_matrix_non_int_N%d_Vg%f_mu%f_T%f_num_freq_pre%d.mat" % (N,Vg,mu,T,num_freq_pre)
	s=foldername + filename
	A = sio.loadmat(s)
	ldos_integrated_non_int_T[ind_T,:,:]=A['ldos_integrated']
	density_non_int_T[ind_T,:,:]=A['density']
	eig_density_non_int_T[ind_T,:]=A['eig_density'][0,:]
	w, v = linalg.eigh(density_non_int_T[ind_T,:,:])
	associated_energies_non_int_T[ind_T,:] = diag(matmul(matmul(v.transpose(),H0),v))
	#Fermi function as reference:

	for ind_f in range(len(wf_fermi)):
	 	fermi_val[ind_T,ind_f] = fermi(wf_fermi[ind_f],0.0,T_list[ind_T])
	
sites=arange(-N,N+1,1)



#Plot stuff:
fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (12,8))

#axis[0].plot(sites,real(diag(H0)),color='black')
#for ind_T in range(len(T_list)):
#	axis[0].plot(sites,real(diag(ldos_integrated_int_T[ind_T,:,:])),linestyle='--')
#axis[0].set_prop_cycle(None)
#for ind_T in range(len(T_list)):
#	axis[0].plot(sites,real(diag(density_int_T[ind_T,:,:])),label=T_list[ind_T])
#legend = axis[0].legend(loc='upper left', fontsize=5)
#plt.gca().set_prop_cycle(None)

for ind_T in range(len(T_list)):
 	#axis[0].plot(arange(0,2*N+1,1),eig_density_int_T[ind_T,:],marker='o', label=T_list[ind_T])
 	#axis[0].plot(associated_energies_int_T[ind_T,:],eig_density_int_T[ind_T,:],marker='o', label=T_list[ind_T],color='blue')
 	axis.plot(associated_energies_int_T_tilde[ind_T,:],eig_density_int_T[ind_T,:],marker='o', label=T_list[ind_T],linestyle='')
legend = axis.legend(loc='upper left', fontsize=5)
plt.gca().set_prop_cycle(None)
for ind_T in range(len(T_list)):
 	#axis[0].plot(arange(0,2*N+1,1),eig_density_non_int_T[ind_T,:],marker='x',linestyle=':')
 	axis.plot(associated_energies_non_int_T[ind_T,:],eig_density_non_int_T[ind_T,:],marker='x',linestyle='')
plt.gca().set_prop_cycle(None)
for ind_T in range(len(T_list)):
 	axis.plot(wf_fermi,fermi_val[ind_T,:])

axis.set_xlim([-3,3])
axis.set_ylim([-0.1,1.1])

savename='density_matrix_int_with_leads_and_full_energies_N%d.pdf' % (2*N+1)
plt.savefig(savename, format = 'pdf')
plt.show();
