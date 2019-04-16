import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex
import Ex_freq_str as Ex_freq_str




#Load normal Flow-Data:
s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation/"
s_file="X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu-1.5000_T0.005000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"

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

h = np.real(ctn.load_as_np(s,"h"))[0,0]

U0 = np.real(ctn.load_as_np(s,"U1"))[0,0]
U1 = np.real(ctn.load_as_np(s,"U2"))[0,0]
T = np.real(ctn.load_as_np(s,"T"))[0,0]

print("U0=",U0)
print("U1=",U1)
print("T=",T)

#[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")


