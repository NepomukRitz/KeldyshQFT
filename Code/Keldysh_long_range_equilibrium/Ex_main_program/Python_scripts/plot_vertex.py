import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex




#Load Flow-Data:
s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
s_file="X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL10_Vg0.2500_h0.000000_mu-1.4700_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no8.mat"
s=s_folder + s_file
L = np.real(ctn.load_as_np(s,"L"))[0,0]
L = L.astype(int)
N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
NfbP = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
NfbP = NfbP.astype(int)
NfbX = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
NfbX = NfbX.astype(int)
wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
Lp_structure = Lp_structure.astype(int)
Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
Lx_structure = Lx_structure.astype(int)
pos_NfbP_2mu = np.real(ctn.load_as_np(s,"pos_NfbP_2mu"))[0,0]
pos_NfbP_2mu = pos_NfbP_2mu.astype(int)

print("pos_NfbP_2mu=",pos_NfbP_2mu)
print("wbP[pos_NfbP_2mu]=",wbP[pos_NfbP_2mu])

[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

#Analyse P structure:

l=1
k=1
ap_r = np.empty([NfbP]) 
ap_i = np.empty([NfbP]) 
for i in range(NfbP):
	Li = Lp_structure[i]
	ap_r[i] = -1.e-7
	ap_i[i] = -1.e-7
	if(abs(l)<=Li and abs(k)<=Li):
		ap_r[i] = np.amax(np.absolute(np.real(aPud_dyn[i][l+Li,k+Li]))) 
		ap_i[i] = np.amax(np.absolute(np.imag(aPud_dyn[i][l+Li,k+Li]))) 
ap_r[pos_NfbP_2mu] = np.amax(np.absolute(np.real(aPud_stat[l+L,k+L]))) 
ap_i[pos_NfbP_2mu] = 0.0 

	
	

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(wbP,ap_r,color='blue',marker='o')
#axis.plot(wbP,ap_i,color='red',marker='o')

axis.set_xlim([-5,5])

plt.show()


