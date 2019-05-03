clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_L5_N10_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.030000_Lambda_ini_10000.000000_Lambda_fin_0.000000.mat')

% Check symmetries of feedback:
diff_feedback=0;
% X-channel:
ind_channel = 4;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - transpose(gamma_data_long_str(ind_channel).m(k_ind,l_ind).m);
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_X_feedback_sym = diff
diff_feedback = max(diff_feedback,diff);
 	
% P-channel equal spin:
ind_channel = 1
diff=0.0;
tmp = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
diff=max(diff,max(max(abs(tmp)))); 

diff_P_feedback_eq_0 = diff
diff_feedback = max(diff_feedback,diff);

ind_channel = 1;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - transpose(gamma_data_long_str(ind_channel).m(k_ind,l_ind).m);
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_P_feedback_eq_sym = diff
diff_feedback = max(diff_feedback,diff);

diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m + gamma_data_long_str(ind_channel).m(2*L+2-l_ind,k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
	
diff_P_feedback_eq_ex_1 = diff
diff_feedback = max(diff_feedback,diff);

diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m + gamma_data_long_str(ind_channel).m(l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
	
diff_P_feedback_eq_ex_2 = diff
diff_feedback = max(diff_feedback,diff);

diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel).m(2*L+2-l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
	
diff_P_feedback_eq_ex_comb = diff
diff_feedback = max(diff_feedback,diff);

diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel+1).m(l_ind,k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_P_feedback_eq_spin = diff
diff_feedback = max(diff_feedback,diff);

%P-channel different spin:
ind_channel = 3;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - transpose(gamma_data_long_str(ind_channel).m(k_ind,l_ind).m);
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_P_feedback_diff_sym = diff
diff_feedback = max(diff_feedback,diff);

%D-channel equal spin:
ind_channel = 5;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - transpose(gamma_data_long_str(ind_channel).m(k_ind,l_ind).m);
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_D_feedback_eq_sym = diff
diff_feedback = max(diff_feedback,diff);

ind_channel = 5;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel).m(2*L+2-l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_D_feedback_eq_ex = diff
diff_feedback = max(diff_feedback,diff);

diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel+1).m(l_ind,k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_D_feedback_eq_spin = diff
diff_feedback = max(diff_feedback,diff);

%D-channel different spin:

ind_channel = 7;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel).m(2*L+2-l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end

diff_D_feedback_diff_ex = diff
diff_feedback = max(diff_feedback,diff);

%Check special B=0 symmetries:

% P-channel
ind_channel = 3;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel).m(2*L+2-l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
diff_P_feedback_diff_ex = diff
diff_feedback = max(diff_feedback,diff);


% X-channel
ind_channel = 4;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - gamma_data_long_str(ind_channel).m(2*L+2-l_ind,2*L+2-k_ind).m;
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
diff_X_feedback_diff_ex = diff
diff_feedback = max(diff_feedback,diff);

% D-channel
ind_channel = 7;
diff=0.0;
for l_ind=1:2*L+1
	for k_ind=1:2*L+1
	 	tmp = gamma_data_long_str(ind_channel).m(l_ind,k_ind).m - transpose(gamma_data_long_str(ind_channel).m(k_ind,l_ind).m);
	 	diff=max(diff,max(max(abs(tmp)))); 
	end
end
diff_D_feedback_diff_sym = diff
diff_feedback = max(diff_feedback,diff);

diff_feedback_gesamt = diff_feedback


