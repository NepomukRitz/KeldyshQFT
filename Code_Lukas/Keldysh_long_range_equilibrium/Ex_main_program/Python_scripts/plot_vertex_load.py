import numpy as np
import cbin_to_numpy as ctn
import matplotlib.pyplot as plt
import Ex_Vertex as Ex_Vertex

def load_complete(s_folder,s_file):
	#Load Flow-Data:
	s=s_folder + s_file
	L = np.real(ctn.load_as_np(s,"L"))[0,0]
	L = L.astype(int)
	N = np.real(ctn.load_as_np(s,"N"))[0,0]
	N = N.astype(int)
	Nff = np.real(ctn.load_as_np(s,"Nff"))[0,0]
	Nff = Nff.astype(int)
	NfbP = np.real(ctn.load_as_np(s,"NfbP"))[0,0]
	NfbP = NfbP.astype(int)
	NfbX = np.real(ctn.load_as_np(s,"NfbX"))[0,0]
	NfbX = NfbX.astype(int)
	wf = np.real(ctn.load_as_np(s,"wf"))[0,:]
	wbP = np.real(ctn.load_as_np(s,"wbP"))[0,:]
	wbX = np.real(ctn.load_as_np(s,"wbX"))[0,:]
	Lp_structure = np.real(ctn.load_as_np(s,"Lp_structure"))[0,:]
	Lp_structure = Lp_structure.astype(int)
	Lx_structure = np.real(ctn.load_as_np(s,"Lx_structure"))[0,:]
	Lx_structure = Lx_structure.astype(int)
	pos_Nff_mu = np.real(ctn.load_as_np(s,"pos_Nff_mu"))[0,0]
	pos_Nff_mu = pos_Nff_mu.astype(int)
	pos_NfbP_2mu = np.real(ctn.load_as_np(s,"pos_NfbP_2mu"))[0,0]
	pos_NfbP_2mu = pos_NfbP_2mu.astype(int)
	pos_NfbX_0 = np.real(ctn.load_as_np(s,"pos_NfbX_0"))[0,0]
	pos_NfbX_0 = pos_NfbX_0.astype(int)
	
	print("pos_NfbP_2mu=",pos_NfbP_2mu)
	print("wbP[pos_NfbP_2mu]=",wbP[pos_NfbP_2mu])
	
	[ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat] = Ex_Vertex.load(s,"gamma")
	return [L,N,Nff,NfbP,NfbX,wf,wbP,wbX,Lp_structure,Lx_structure,pos_Nff_mu,pos_NfbP_2mu,pos_NfbX_0,ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat]

#s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
#s_file="X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL2_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no4.mat"
#[L0,N0,Nff0,NfbP0,NfbX0,wf0,wbP0,wbX0,Lp_structure0,Lx_structure0,pos_Nff_mu0,pos_NfbP_2mu0,pos_NfbX_00,ERetu0, ERetd0, aPuu_dyn0, aPuu_stat0, aPdd_dyn0, aPdd_stat0, aPud_dyn0, aPud_stat0, aXud_dyn0, aXud_stat0, aDuu_dyn0, aDuu_stat0, aDdd_dyn0, aDdd_stat0, aDud_dyn0, aDud_stat0] = load_complete(s_folder,s_file)
#
#s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
#s_file="X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL5_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no10.mat"
#[L1,N1,Nff1,NfbP1,NfbX1,wf1,wbP1,wbX1,Lp_structure1,Lx_structure1,pos_Nff_mu1,pos_NfbP_2mu1,pos_NfbX_01,ERetu1, ERetd1, aPuu_dyn1, aPuu_stat1, aPdd_dyn1, aPdd_stat1, aPud_dyn1, aPud_stat1, aXud_dyn1, aXud_stat1, aDuu_dyn1, aDuu_stat1, aDdd_dyn1, aDdd_stat1, aDud_dyn1, aDud_stat1] = load_complete(s_folder,s_file)
#
#s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
#s_file="X_L5_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL10_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no8.mat"
#[L2,N2,Nff2,NfbP2,NfbX2,wf2,wbP2,wbX2,Lp_structure2,Lx_structure2,pos_Nff_mu2,pos_NfbP_2mu2,pos_NfbX_02,ERetu2, ERetd2, aPuu_dyn2, aPuu_stat2, aPdd_dyn2, aPdd_stat2, aPud_dyn2, aPud_stat2, aXud_dyn2, aXud_stat2, aDuu_dyn2, aDuu_stat2, aDdd_dyn2, aDdd_stat2, aDud_dyn2, aDud_stat2] = load_complete(s_folder,s_file)



s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
s_file="X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
[L0,N0,Nff0,NfbP0,NfbX0,wf0,wbP0,wbX0,Lp_structure0,Lx_structure0,pos_Nff_mu0,pos_NfbP_2mu0,pos_NfbX_00,ERetu0, ERetd0, aPuu_dyn0, aPuu_stat0, aPdd_dyn0, aPdd_stat0, aPud_dyn0, aPud_stat0, aXud_dyn0, aXud_stat0, aDuu_dyn0, aDuu_stat0, aDdd_dyn0, aDdd_stat0, aDud_dyn0, aDud_stat0] = load_complete(s_folder,s_file)

s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
s_file="X_L0_Lu0_N30_Nff1500_NfbP500_NfbX500_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
[L1,N1,Nff1,NfbP1,NfbX1,wf1,wbP1,wbX1,Lp_structure1,Lx_structure1,pos_Nff_mu1,pos_NfbP_2mu1,pos_NfbX_01,ERetu1, ERetd1, aPuu_dyn1, aPuu_stat1, aPdd_dyn1, aPdd_stat1, aPud_dyn1, aPud_stat1, aXud_dyn1, aXud_stat1, aDuu_dyn1, aDuu_stat1, aDdd_dyn1, aDdd_stat1, aDud_dyn1, aDud_stat1] = load_complete(s_folder,s_file)

s_folder="/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz/"
s_file="X_L0_Lu0_N30_Nff1500_NfbP250_NfbX250_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4600_T0.000000_Uc0.50_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
[L2,N2,Nff2,NfbP2,NfbX2,wf2,wbP2,wbX2,Lp_structure2,Lx_structure2,pos_Nff_mu2,pos_NfbP_2mu2,pos_NfbX_02,ERetu2, ERetd2, aPuu_dyn2, aPuu_stat2, aPdd_dyn2, aPdd_stat2, aPud_dyn2, aPud_stat2, aXud_dyn2, aXud_stat2, aDuu_dyn2, aDuu_stat2, aDdd_dyn2, aDdd_stat2, aDud_dyn2, aDud_stat2] = load_complete(s_folder,s_file)
