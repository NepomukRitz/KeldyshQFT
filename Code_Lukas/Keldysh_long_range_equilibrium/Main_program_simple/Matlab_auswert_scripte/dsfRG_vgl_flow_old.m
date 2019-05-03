clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_vgl_flow_old.mat')

for ind=1:length(y_old(3).m)
 	%A(ind) = y_old(5).m(ind).m(11,11);
 	A(ind) = y_old(3).m(ind).m(11,11) - U(11,11)/4;
end

for ind=1:length(gamma_data_short_str(3).m)
 	B(ind) = gamma_data_short_str(3).m(ind).m(11,11);
end

max(max(abs(A-B)))

figure
hold 
plot(wf,imag(A))
plot(wf,imag(B))
xlim([-10 10])

figure
hold 
plot(wf,real(A))
plot(wf,real(B))
xlim([-10 10])

%Check Consistency between long and short structure:

%for ind_channel=1:3
%	P_short = gamma_data_short_str(ind_channel).m(pos_NfbP_2mu+1).m;
%	P_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
%	diff_consistency(ind_channel) = max(max(abs(P_short - P_long)));
%end
for ind_channel=4:4
 	X_short = gamma_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
	X_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
	diff_consistency(ind_channel) = max(max(abs(X_short - X_long)));
end

%figure
%imagesc(real(gamma_data_short_str(4).m(pos_NfbX_0+1).m) - gamma_data_long_str(4).m(L+1,L+1).m)
%debug_real= max(max(abs(  real(gamma_data_short_str(4).m(pos_NfbX_0+1).m)  )))
%
%for ind_channel=5:5
% 	D_short = gamma_data_short_str(ind_channel).m(pos_NfbX_0+1).m;
%	D_long  = gamma_data_long_str(ind_channel).m(L+1,L+1).m;
%	diff_consistency(ind_channel) = max(max(abs(D_short - D_long)));
%end

diff_complete = max(diff_consistency)

%Debug part

%debug_real= max(max(abs(  real(gamma_data_short_str(4).m(pos_NfbX_0+1).m)  )))
figure
imagesc(real(gamma_data_short_str(4).m(pos_NfbX_0+1).m))
figure
imagesc(real(gamma_data_long_str(4).m(L+1,L+1).m))
figure
imagesc(real(gamma_data_long_str(4).m(L+1,L+1).m./gamma_data_short_str(4).m(pos_NfbX_0+1).m))


%	%Check at -infinity
%	A_minf = y_old(3).m(1).m - U/4
%	B_minf = gamma_data_short_str(3).m(1).m;
%	
%	A_pinf = y_old(3).m(end).m - U/4
%	B_pinf = gamma_data_short_str(3).m(end).m;
%	
%Check all elements:

for ind=1:length(y_old(3).m)
 	A = y_old(6).m(ind).m;
 	%A = y_old(4).m(ind).m - U/4;
 	B = gamma_data_short_str(6).m(ind).m;
	diff(ind) = max(max(abs(A-B )));
	%diff(ind) = max(max(abs(A-B - (A_minf - B_minf) )));
end
diff_max = max(diff)

figure
plot(wf, diff)
%xlim([-10 10])




%Check Verschwinden von Puu, Pdd, Dud:
%
%for ind=1:length(gamma_data_short_str(3).m)
% 	Puu_short_max(ind) = max(max( abs( gamma_data_short_str(1).m(ind).m ) ) ); 
% 	Pdd_short_max(ind) = max(max( abs( gamma_data_short_str(2).m(ind).m ) ) );
% 	Dud_short_max(ind) = max(max( abs( gamma_data_short_str(7).m(ind).m ) ) );
%end
%
%Puu_short_ges = max(Puu_short_max)
%Pdd_short_ges = max(Pdd_short_max)
%Dud_short_ges = max(Dud_short_max)
%
%for ind1=1:length(gamma_data_long_str(3).m)
%	for ind2=1:length(gamma_data_long_str(3).m)
% 		Puu_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(1).m(ind1,ind2).m ) ) ); 
% 		Pdd_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(2).m(ind1,ind2).m ) ) );
% 		Dud_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(7).m(ind1,ind2).m ) ) );
%	end
%
%end
%
%Puu_long_ges = max(Puu_long_max)
%Pdd_long_ges = max(Pdd_long_max)
%Dud_long_ges = max(Dud_long_max)
