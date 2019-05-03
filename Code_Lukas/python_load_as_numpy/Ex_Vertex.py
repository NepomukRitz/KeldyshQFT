import numpy as np 
import cbin_to_numpy as ctn
import Ex_freq_str as Ex_freq_str

def load(filename,variable):
	s_ERetu = variable + "_ERetu"
	s_ERetd = variable + "_ERetd"
	s_aPuu = variable + "_aPuu"
	s_aPdd = variable + "_aPdd"
	s_aPud = variable + "_aPud"
	s_aXud = variable + "_aXud"
	s_aDuu = variable + "_aDuu"
	s_aDdd = variable + "_aDdd"
	s_aDud = variable + "_aDud"
	ERetu = ctn.load_as_np(filename,s_ERetu)
	ERetd = ctn.load_as_np(filename,s_ERetd)
	[L, N, L_structure, wb, aPuu_dyn, aPuu_stat] = Ex_freq_str.load(filename,s_aPuu)
	[L, N, L_structure, wb, aPdd_dyn, aPdd_stat] = Ex_freq_str.load(filename,s_aPdd)
	[L, N, L_structure, wb, aPud_dyn, aPud_stat] = Ex_freq_str.load(filename,s_aPud)
	[L, N, L_structure, wb, aXud_dyn, aXud_stat] = Ex_freq_str.load(filename,s_aXud)
	[L, N, L_structure, wb, aDuu_dyn, aDuu_stat] = Ex_freq_str.load(filename,s_aDuu)
	[L, N, L_structure, wb, aDdd_dyn, aDdd_stat] = Ex_freq_str.load(filename,s_aDdd)
	[L, N, L_structure, wb, aDud_dyn, aDud_stat] = Ex_freq_str.load(filename,s_aDud)
	return [ERetu, ERetd, aPuu_dyn, aPuu_stat, aPdd_dyn, aPdd_stat, aPud_dyn, aPud_stat, aXud_dyn, aXud_stat, aDuu_dyn, aDuu_stat, aDdd_dyn, aDdd_stat, aDud_dyn, aDud_stat]


