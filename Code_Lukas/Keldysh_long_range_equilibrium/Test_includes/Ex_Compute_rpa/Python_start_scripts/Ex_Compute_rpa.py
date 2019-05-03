import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex


#Load RPA-Data:	
s_folder="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Compute_rpa/"
s_file="Ex_Compute_rpa.mat"

s=s_folder + s_file

wbP_rpa = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX_rpa = np.real(ctn.load_as_np(s,"wbX"))[0,:]
[ERetu_rpa, ERetd_rpa, aPuu_dyn_rpa, aPuu_stat_rpa, aPdd_dyn_rpa, aPdd_stat_rpa, aPud_dyn_rpa, aPud_stat_rpa, aXud_dyn_rpa, aXud_stat_rpa, aDuu_dyn_rpa, aDuu_stat_rpa, aDdd_dyn_rpa, aDdd_stat_rpa, aDud_dyn_rpa, aDud_stat_rpa] = Ex_Vertex.load(s,"gamma_rpa")



#Load Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG/"
s_file="ex_dsfRG_rpa_L2_Lu2_N5_Nff_10_NfbP_10_NfbX_10_num_freq_pre_30000_Vg_0.250000_h_0.100000_mu_-1.475000_T_0.010000_U0_0.300000_U1_0.100000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_1.mat"

s=s_folder + s_file

wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

print("diff_puu=",np.amax(np.absolute(aPuu_dyn_rpa - aPuu_dyn)))
print("diff_pdd=",np.amax(np.absolute(aPdd_dyn_rpa - aPdd_dyn)))
print("diff_pud=",np.amax(np.absolute(aPud_dyn_rpa - aPud_dyn)))
print("diff_xud=",np.amax(np.absolute(aXud_dyn_rpa - aXud_dyn)))
print("diff_duu=",np.amax(np.absolute(aDuu_dyn_rpa - aDuu_dyn)))
print("diff_ddd=",np.amax(np.absolute(aDdd_dyn_rpa - aDdd_dyn)))
print("diff_dud=",np.amax(np.absolute(aDud_dyn_rpa - aDud_dyn)))

print("abs_puu=",np.amax(np.absolute(aPuu_dyn)))
print("abs_puu_rpa=",np.amax(np.absolute(aPuu_dyn_rpa)))

print("abs_pdd=",np.amax(np.absolute(aPdd_dyn)))
print("abs_pdd_rpa=",np.amax(np.absolute(aPdd_dyn_rpa)))

print("abs_pud=",np.amax(np.absolute(aPud_dyn)))
print("abs_pud_rpa=",np.amax(np.absolute(aPud_dyn_rpa)))

print("abs_xud=",np.amax(np.absolute(aXud_dyn)))
print("abs_xud_rpa=",np.amax(np.absolute(aXud_dyn_rpa)))

print("abs_duu=",np.amax(np.absolute(aDuu_dyn)))
print("abs_duu_rpa=",np.amax(np.absolute(aDuu_dyn_rpa)))

print("abs_ddd=",np.amax(np.absolute(aDdd_dyn)))
print("abs_ddd_rpa=",np.amax(np.absolute(aDdd_dyn_rpa)))

print("abs_dud=",np.amax(np.absolute(aDud_dyn)))
print("abs_dud_rpa=",np.amax(np.absolute(aDud_dyn_rpa)))


aPuu_dyn_rpa_abs = np.empty([np.size(wbP)])
aPuu_dyn_abs = np.empty([np.size(wbP)])
for i in range(np.size(wbP)):
	aPuu_dyn_rpa_abs[i] = np.amax(np.absolute(aPuu_dyn_rpa[i,:,:,:,:]))
	aPuu_dyn_abs[i] = np.amax(np.absolute(aPuu_dyn[i,:,:,:,:]))

fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
axis.plot(wbX_rpa,np.real(aDuu_dyn_rpa[:,0,0,1,1]),color='blue',marker='o')
axis.plot(wbX_rpa,np.imag(aDuu_dyn_rpa[:,0,0,1,1]),color='red',marker='o')
axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='blue',linestyle='--',marker='o')
axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='red',linestyle='--',marker='o')



#	
#axis.plot(wbX_rpa,np.real(aXud_dyn_rpa[:,0,0,2,2]),color='blue',marker='o')
#axis.plot(wbX_rpa,np.imag(aXud_dyn_rpa[:,0,0,2,2]),color='red',marker='o')
#axis.plot(wbX,np.real(aXud_dyn[:,0,0,2,2]),color='black',linestyle='--',marker='o')
#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,2,2]),color='black',linestyle='--',marker='o')
#	
axis.set_xlim([-6,6])

plt.show()


