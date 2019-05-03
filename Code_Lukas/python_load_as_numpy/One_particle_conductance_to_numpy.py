import numpy as np 
import cbin_to_numpy as ctn
import One_particle_conductance_for_py as Opcfp 

def compute(spin,mu,h,self_up,self_down,Lambda):
	cond_obj = Opcfp.Conductance_for_py_py()
	[Nges, dump] = np.shape(self_up)
	cond_obj.set_parameters(Nges,mu,h,Lambda)
	for i in range(Nges):
		for j in range(Nges):
			cond_obj.populate_ERet_at_mu(i,j,self_up[i,j],self_down[i,j])
	return cond_obj.compute(spin)

def get_cond_curve_non_int(spin,mu_list,h,H):
	Lambda=1e-8
	N_mu = np.size(mu_list)
	result = np.empty([N_mu])
	for i in range(N_mu):
	 	mu = mu_list[i]
		result[i] = compute(spin,mu,h,H,H,Lambda)
	return result

def cond_combined(s,mu_list):
	H = np.real(ctn.load_as_np(s,"hamiltonian"))
	h = np.real(ctn.load_as_np(s,"h"))[0,0]
	cond_up_list = get_cond_curve_non_int(1,mu_list,h,H)
	cond_down_list = get_cond_curve_non_int(0,mu_list,h,H)
	return cond_up_list + cond_down_list

