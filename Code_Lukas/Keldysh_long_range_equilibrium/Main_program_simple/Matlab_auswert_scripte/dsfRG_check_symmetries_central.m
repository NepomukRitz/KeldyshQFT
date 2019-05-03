clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_L5_N10_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.030000_Lambda_ini_10000.000000_Lambda_fin_0.000000.mat')

% Check symmetries of central part:

%P-channel:

ind_channel = 1
diff=0.0;
for ind_f=1:NfbP
	tmp = gamma_data_short_str(ind_channel).m(ind_f).m;
	diff=max(diff,max(max(abs(tmp)))); 
end

diff_P_central_eq_0 = diff

ind_channel = 3
diff=0.0;
for ind_f=1:NfbP
	tmp = gamma_data_short_str(ind_channel).m(ind_f).m - transpose(gamma_data_short_str(ind_channel).m(ind_f).m);
	diff=max(diff,max(max(abs(tmp)))); 
end

diff_P_central_diff_sym = diff

%X-channel:

ind_channel = 4
diff=0.0;
for ind_f=1:NfbX
	tmp = gamma_data_short_str(ind_channel).m(ind_f).m - transpose(gamma_data_short_str(ind_channel).m(ind_f).m);
	diff=max(diff,max(max(abs(tmp)))); 
end

diff_X_central_diff_sym = diff

diff=0.0;
for ind_j=1:2*N+1
	for ind_i=1:2*N+1
	 	for ind_f=1:NfbX
		 	tmp_1(ind_f) = conj(gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i));
			tmp_2(ind_f) = gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i);
		end
		A_tmp = interp1(wbX, tmp_1,-wbX);
		tmp = A_tmp - tmp_2;
		diff=max(diff,max(abs(tmp)));
	end
end

diff_X_central_diff_ex = diff

%D-channel:

ind_channel = 5
diff=0.0;
for ind_f=1:NfbX
	tmp = gamma_data_short_str(ind_channel).m(ind_f).m - transpose(gamma_data_short_str(ind_channel).m(ind_f).m);
	diff=max(diff,max(max(abs(tmp)))); 
end

diff_D_central_eq_sym = diff

diff=0.0;
for ind_j=1:2*N+1
	for ind_i=1:2*N+1
	 	for ind_f=1:NfbX
		 	tmp_1(ind_f) = conj(gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i));
			tmp_2(ind_f) = gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i);
		end
		A_tmp = interp1(wbX, tmp_1,-wbX);
		tmp = A_tmp - tmp_2;
		diff=max(diff,max(abs(tmp)));
	end
end

diff_D_central_eq_ex = diff

ind_channel = 7
diff=0.0;
for ind_j=1:2*N+1
	for ind_i=1:2*N+1
	 	for ind_f=1:NfbX
		 	tmp_1(ind_f) = conj(gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i));
			tmp_2(ind_f) = gamma_data_short_str(ind_channel).m(ind_f).m(ind_j,ind_i);
		end
		A_tmp = interp1(wbX, tmp_1,-wbX);
		tmp = A_tmp - tmp_2;
		diff=max(diff,max(abs(tmp)));
	end
end

diff_D_central_diff_ex = diff


