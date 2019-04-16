import numpy as np 
import cbin_to_numpy as ctn

def load(filename,variable):
	s_L =  variable + "_L"
	s_N =  variable + "_N"
	s_L_structure =  variable + "_L_structure"
	s_wb =  variable + "_wb"
	s_dyn = variable + "_dynamic_str"
	s_stat = variable + "_static_str"

	tmp = ctn.load_as_np(filename,s_L) 
	L=int(np.real(tmp[0,0]))
	tmp = ctn.load_as_np(filename,s_N) 
	N=int(np.real(tmp[0,0]))
	tmp = ctn.load_as_np(filename,s_L_structure) 
	L_structure = np.real(tmp[0,:])
	L_structure = L_structure.astype(int)
	tmp = ctn.load_as_np(filename,s_wb) 
	wb = np.real(tmp[0,:])
	tmp = ctn.load_as_np(filename,s_dyn) 
	dyn_raw = tmp[0,:]
	tmp = ctn.load_as_np(filename,s_stat) 
	stat_raw = np.real(tmp[0,:])

	N_freq = np.size(L_structure)
	Nges = 2*N+1
	Lges = 2*L+1
	dyn=[]
	z=0
	for i in range(N_freq):
		Li = L_structure[i]
		Liges = 2*Li+1
		list_l=[]
		for l in range(Liges):
			list_k =[]
			for k in range(Liges):
				tmp = np.empty([Nges-abs(l-Li),Nges-abs(k-Li)], dtype = np.complex128)
				for j1 in range(Nges-abs(l-Li)):
					for j2 in range(Nges-abs(k-Li)):
						tmp[j1,j2] = dyn_raw[z]
						z=z+1
				list_k.append(tmp)
			list_l.append(list_k)
		tmp_l = np.asarray(list_l)
		dyn.append(tmp_l)
	dyn = np.asarray(dyn)

	z=0
	list_l = []
	for l in range(Lges):
		list_k = []
		for k in range(Lges):
			tmp = np.empty([Nges-abs(l-L),Nges-abs(k-L)], dtype = np.float64)
			for j1 in range(Nges-abs(l-L)):
				for j2 in range(Nges-abs(k-L)):
					tmp[j1,j2] = stat_raw[z]
					z=z+1
			list_k.append(tmp)
		list_l.append(list_k)
	stat = np.asarray(list_l)
	return [L, N, L_structure, wb, dyn, stat]

def abs_dyn(dyn,L_structure):
	ret = 0.0
	N_freq=np.size(L_structure)
	for i in range(N_freq):
		for l in range(2*L_structure[i]+1):
			for k in range(2*L_structure[i]+1):
				ret = max(ret,np.amax(np.absolute(dyn[i][l,k])))
	return ret

def abs_stat(stat,L):
	ret = 0.0
	for l in range(2*L+1):
		for k in range(2*L+1):
			ret = max(ret,np.amax(np.absolute(stat[l,k])))
	return ret
		
