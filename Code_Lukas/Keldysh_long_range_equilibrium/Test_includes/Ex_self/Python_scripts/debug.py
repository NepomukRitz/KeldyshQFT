import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex




#Load debug :
s_folder="../Start_scripts/"
s_file="debug.mat"


s=s_folder + s_file

N = np.real(ctn.load_as_np(s,"N"))[0,0]
N = N.astype(int)
#wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
#wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
#wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
#[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")

Su = ctn.load_as_np(s,"Su")
freq_pre = np.real(ctn.load_as_np(s,"freq_pre")[0,:])



fig, axis = plt.subplots(nrows = 1, ncols = 1, figsize = (6,3))
##	#axis.plot(wbX_old,np.real(aDuu_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbX_old,np.imag(aDuu_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbX,np.real(aDuu_dyn[:,0,0,1,1]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbX,np.imag(aDuu_dyn[:,0,0,1,1]),color='orange',linestyle='--',marker='o')
##	
##	#axis.plot(wbP_old,np.real(aPud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbP_old,np.imag(aPud_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbP,np.real(aPud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbP,np.imag(aPud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
##	
##	#axis.plot(wbX_old,np.real(aXud_dyn_old[:,0,0,N,N]),color='blue',marker='o')
##	#axis.plot(wbX_old,np.imag(aXud_dyn_old[:,0,0,N,N]),color='red',marker='o')
##	#axis.plot(wbX,np.real(aXud_dyn[:,0,0,N,N]),color='black',linestyle='--',marker='o')
##	#axis.plot(wbX,np.imag(aXud_dyn[:,0,0,N,N]),color='orange',linestyle='--',marker='o')
##	
print(Su[0,0,N,N])
axis.plot(freq_pre,np.absolute(Su[0,:,N,N]),color='blue')
axis.plot(freq_pre,np.absolute(Su[0,:,0,0]),color='red')
#axis.plot(wf_old,np.imag(ERetu_old[0,:,N,N]),color='red',marker='o')
#axis.plot(wf,np.real(ERetu[0,:,N,N]),color='black',linestyle='--',marker='o')
#axis.plot(wf,np.imag(ERetu[0,:,N,N]),color='orange',linestyle='--',marker='o')
#
#axis.set_xlim([-5,5])
#
plt.show()


