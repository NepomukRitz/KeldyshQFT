import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex

def absolute(A,Nff,N,L_structure):
	value=np.empty([Nff])
	for i in range(Nff):
	 	L_inner = L_structure[i]
		value[i] = np.amax(np.absolute(A[i][L_inner,L_inner][N,N]))
	return value
	 	

#Load RPA-Data:	
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Compute_rpa/"
s_file="Ex_Compute_rpa.mat"

s=s_folder + s_file

wbP_rpa = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
[ERetu_rpa, ERetd_rpa, aPuu_dyn_rpa, aPuu_stat_rpa, aPdd_dyn_rpa, aPdd_stat_rpa, aPud_dyn_rpa, aPud_stat_rpa, aXud_dyn_rpa, aXud_stat_rpa, aDuu_dyn_rpa, aDuu_stat_rpa, aDdd_dyn_rpa, aDdd_stat_rpa, aDud_dyn_rpa, aDud_stat_rpa] = Ex_Vertex.load(s,"gamma_rpa")

N_rpa = np.real(ctn.load_as_np(s,"N"))[0,0]
N_rpa = N_rpa.astype(int)
L_rpa = np.real(ctn.load_as_np(s,"L"))[0,0]
L_rpa = L_rpa.astype(int)
wbP_rpa = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf_rpa = np.real(ctn.load_as_np(s,"wf"))[0,:]
Lp_structure_rpa = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
Lp_structure_rpa = Lp_structure_rpa.astype(int)
Lx_structure_rpa = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
Lx_structure_rpa = Lx_structure_rpa.astype(int)
NfbP_rpa = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
NfbP_rpa = NfbP_rpa.astype(int)
NfbX_rpa = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
NfbX_rpa = NfbX_rpa.astype(int)



#Load Flow-Data:
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_dsfRG/"
s_file="utof_debug.mat"

s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
L = np.real(ctn.load_as_np(s,"L"))[0,0]
L = L.astype(int)
wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
Lp_structure = Lp_structure.astype(int)
Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
Lx_structure = Lx_structure.astype(int)
NfbP = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
NfbP = NfbP.astype(int)
NfbX = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
NfbX = NfbX.astype(int)

print(Lp_structure)
print(Lx_structure)

[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

diff = 0
for i in range(NfbP):
	L = Lp_structure[i]
	for l in range(-L,L+1):
		for k in range(-L,L+1):
			tmp = aPuu_dyn[i][l+L,k+L]- aPuu_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
			tmp = aPdd_dyn[i][l+L,k+L]- aPdd_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
			tmp = aPud_dyn[i][l+L,k+L]- aPud_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))

print("diff P=",diff)

diff = 0
for i in range(NfbX):
	L = Lx_structure[i]
	for l in range(-L,L+1):
		for k in range(-L,L+1):
			tmp = aXud_dyn[i][l+L,k+L]- aXud_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
			tmp = aDuu_dyn[i][l+L,k+L]- aDuu_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
			tmp = aDdd_dyn[i][l+L,k+L]- aDdd_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
			tmp = aDud_dyn[i][l+L,k+L]- aDud_dyn_rpa[i][l+L,k+L]
			diff = max(diff, np.amax(np.absolute(tmp)))
print("diff XD=",diff)

fig, axis = plt.subplots(nrows = 1, ncols = 2, figsize = (6,3))

axis[0].plot(wbP,absolute(aPud_dyn,NfbP,N,Lp_structure),color='blue',marker='o')
axis[0].plot(wbP_rpa,absolute(aPud_dyn_rpa,NfbP_rpa,N_rpa,Lp_structure_rpa),color='red',marker='o')

axis[1].plot(wbX,absolute(aXud_dyn,NfbX,N,Lx_structure),color='blue',marker='o')
axis[1].plot(wbX_rpa,absolute(aXud_dyn_rpa,NfbX_rpa,N_rpa,Lx_structure_rpa),color='red',marker='o')

axis[0].set_xlim([-6,6])
axis[1].set_xlim([-6,6])
plt.show()


