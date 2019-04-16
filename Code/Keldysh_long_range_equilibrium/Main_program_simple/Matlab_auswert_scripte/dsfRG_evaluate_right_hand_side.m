clear all
close all


load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_evaluate_right_hand_side.mat')


%Check Consistency between long and short structure:

for ind_channel=1:3
	P_short = dy_new_data_short_str(ind_channel).m(pos_NfbP_2mu+1).m;
	P_long  = dy_new_data_long_str(ind_channel).m(L+1,L+1).m;
	diff_consistency(ind_channel) = max(max(abs(P_short - P_long)));
end
for ind_channel=4:4
 	X_short = dy_new_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
	X_long  = dy_new_data_long_str(ind_channel).m(L+1,L+1).m;
	diff_consistency(ind_channel) = max(max(abs(X_short - X_long)));
end
for ind_channel=5:7
 	D_short = dy_new_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
	D_long  = dy_new_data_long_str(ind_channel).m(L+1,L+1).m;
	diff_consistency(ind_channel) = max(max(abs(D_short - D_long)));
end

diff_complete = max(diff_consistency)


%	for ind_channel=1:3
%		P_short = gamma_data_short_str(ind_channel).m(pos_NfbP_2mu+1).m;
%		P_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
%		diff_consistency(ind_channel) = max(max(abs(P_short - P_long)));
%	end
%	for ind_channel=4:4
%	 	X_short = gamma_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
%		X_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
%		diff_consistency(ind_channel) = max(max(abs(X_short - X_long)));
%	end
%	for ind_channel=5:7
%	 	D_short = gamma_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
%		D_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
%		diff_consistency(ind_channel) = max(max(abs(D_short - D_long)));
%	end
%	
%	diff_complete = max(diff_consistency)



%	%Check magnetic field zero:
%	
%	diff=0;
%	for ind_wbP=1:length(wbP)
%	 	tmp = max(max(abs(  gamma_data_short_str(1).m(ind_wbP).m - gamma_data_short_str(2).m(ind_wbP).m  )));
%		diff = max(diff,tmp);
%	 	tmp = max(max(abs(  gamma_data_short_str(5).m(ind_wbP).m - gamma_data_short_str(6).m(ind_wbP).m  )));
%		diff = max(diff,tmp);
%	end
%	diff_mag_short = diff
%	
%	diff=0;
%	for ind_l=1:2*L+1
%	 	for ind_k=1:2*L+1
%			tmp = max(max(abs(  gamma_data_long_str(1).m(ind_l,ind_k).m - gamma_data_long_str(2).m(ind_l,ind_k).m  )));
%			diff = max(diff,tmp);
%			tmp = max(max(abs(  gamma_data_long_str(5).m(ind_l,ind_k).m - gamma_data_long_str(6).m(ind_l,ind_k).m  )));
%			diff = max(diff,tmp);
%		end
%	end
%	diff_mag_long = diff
%	
%	%Check Particle exchange:
%	
%	diff=0;
%	for ind_l=1:2*L+1
%	 	for ind_k=1:2*L+1
%			tmp = max(max(abs(  gamma_data_long_str(1).m(ind_l,ind_k).m + gamma_data_long_str(1).m(2*L+2-ind_l,ind_k).m  )));
%			diff = max(diff,tmp);
%			tmp = max(max(abs(  gamma_data_long_str(1).m(ind_l,ind_k).m + gamma_data_long_str(1).m(ind_l,2*L+2-ind_k).m  )));
%			diff = max(diff,tmp);
%			tmp = max(max(abs(  gamma_data_long_str(1).m(ind_l,ind_k).m - gamma_data_long_str(1).m(2*L+2-ind_l,2*L+2-ind_k).m  )));
%			diff = max(diff,tmp);
%		end
%	end
%	diff_p_exchange_particles_long = diff
%	
%	%Debug part
%	debug_real= max(max(abs(  real(gamma_data_short_str(ind_channel).m(pos_NfbX_0+1).m)  )))
%	figure
%	imagesc(real(gamma_data_short_str(5).m(pos_NfbX_0+1).m) - gamma_data_long_str(5).m(L+1,L+1).m)


