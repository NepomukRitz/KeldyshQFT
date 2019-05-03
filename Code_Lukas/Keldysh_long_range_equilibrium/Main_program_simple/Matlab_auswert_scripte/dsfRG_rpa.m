clear all
close all
load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_rpa_flow_L0_N10_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.100000_T_0.030000_Lambda_ini_10000.000000_Lambda_fin_0.000000.mat')

for ind=1:length(wf)
 	Sigma(ind) = gamma_data_short_str(8).m(ind).m(N+1,N+1);
	ap_ud(ind) = gamma_data_short_str(3).m(ind).m(N+1,N+1);
	ax_ud(ind) = gamma_data_short_str(4).m(ind).m(N+1,N+1);
	ad_ud(ind) = gamma_data_short_str(7).m(ind).m(N+1,N+1);
end
figure
hold all
plot(wf,Sigma)
plot(wbP,imag(ap_ud))	
plot(wbX,imag(ax_ud))	
plot(wbX,imag(ad_ud))	
xlim([-4 4])
