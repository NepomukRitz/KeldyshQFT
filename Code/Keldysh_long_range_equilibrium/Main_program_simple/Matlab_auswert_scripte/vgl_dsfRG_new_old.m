clear all
close all

New = load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_L0_N10_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.030000_Lambda_ini_10000.000000_Lambda_fin_0.000000.mat')

for ind=1:length(New.gamma_data_short_str(3).m)
 	B(ind) = New.gamma_data_short_str(5).m(ind).m(11,11);
end

Old = load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/flow_old_L0_N10_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.030000_Lambda_ini_10000.000000_Lambda_fin_0.000000.mat')

for ind=1:length(Old.y_old(3).m)
 	A(ind) = Old.y_old(5).m(ind).m(11,11);
 	%A(ind) = y_old(4).m(ind).m(11,11) - U(11,11)/4;
end

max(max(abs(A-B)))

figure
hold 
plot(Old.wf,imag(A))
plot(New.wf,imag(B))
xlim([-10 10])

figure
hold 
plot(Old.wf,real(A))
plot(New.wf,real(B))
xlim([-10 10])

%Check Consistency between long and short structure:

for ind_channel=1:3
	P_short = New.gamma_data_short_str(ind_channel).m(New.pos_NfbP_2mu+1).m;
	P_long  = New.gamma_data_long_str(ind_channel).m(New.L+1,New.L+1).m;
	diff_consistency(ind_channel) = max(max(abs(P_short - P_long)));
end
for ind_channel=4:4
 	X_short = New.gamma_data_short_str(ind_channel).m(New.pos_NfbX_0+1).m;
	X_long  = New.gamma_data_long_str(ind_channel).m(New.L+1,New.L+1).m;
	diff_consistency(ind_channel) = max(max(abs(X_short - X_long)));
end
for ind_channel=5:7
 	D_short = New.gamma_data_short_str(ind_channel).m(New.pos_NfbX_0+1).m;
	D_long  = New.gamma_data_long_str(ind_channel).m(New.L+1,New.L+1).m;
	diff_consistency(ind_channel) = max(max(abs(D_short - D_long)));
end

diff_complete = max(diff_consistency)
