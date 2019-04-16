import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex




#Load Flow-Data:
s_folder="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_Precomputation/"
s_file="preintegrate.mat"

s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)

Su_ret_integrated_slow = ctn.load_as_np(s,"Su_ret_integrated_slow")
Su_ret_integrated_fast = ctn.load_as_np(s,"Su_ret_integrated_fast")
Su_kel_integrated_slow = ctn.load_as_np(s,"Su_kel_integrated_slow")
Su_kel_integrated_fast = ctn.load_as_np(s,"Su_kel_integrated_fast")
Su = ctn.load_as_np(s,"Su")

print(np.amax(np.absolute(Su_ret_integrated_slow - Su_ret_integrated_fast)))



fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
#	#axis.plot(wbX_old,np.real(aDuu_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_old,np.imag(aDuu_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbP_old,np.real(aPud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbP_old,np.imag(aPud_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbP,np.real(aPud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbP,np.imag(aPud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	
#	#axis.plot(wbX_old,np.real(aXud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
#	#axis.plot(wbX_old,np.imag(aXud_dyn_old[:,0,0,N,N]),color='red',marker='o')
#	#axis.plot(wbX,np.real(aXud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
#	#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
#	

#axis.plot(range(0,30000),np.imag(Su_ret_integrated_slow[0,:,N,N]),color='green',marker='x',linestyle='--')
#axis.plot(range(0,30000),np.imag(Su_ret_integrated_fast[0,:,N,N]),color='red',marker='x',linestyle='--')

#axis.plot(range(0,30000),np.imag(Su_kel_integrated_fast[0,:,N,N] - Su_kel_integrated_slow[0,:,N,N] ),color='red',marker='x',linestyle='--')
axis.plot(range(0,30000),np.imag(Su_kel_integrated_slow[0,:,N,N]),color='red',marker='x',linestyle='--')

#axis.set_xlim([0,10])

plt.show()


