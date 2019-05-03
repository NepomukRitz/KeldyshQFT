clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_vgl_right_hand_side.mat')

for ind=1:length(wf)
 	A_old(ind) = dy_old(6).m(ind).m(31,31);
 	A_new(ind) = dy_new_data_short_str(6).m(ind).m(31,31);
end

figure
hold 
plot(wf,imag(A_old))
plot(wf,imag(A_new))
plot(wf,real(A_old))
plot(wf,real(A_new))

xlim([-10 10])

max(abs(A_old - A_new))

% Kompletter check:

maximum = 0.0;
for ind=1:length(wf)
 	A_old = dy_old(6).m(ind).m;
 	A_new = dy_new_data_short_str(6).m(ind).m;
	maximum = max(maximum,max(max( abs(A_old - A_new))));
	diff(ind) = max(max( abs(A_old - A_new)));
	%max_vec_real(ind) = max(max(real(A_old - A_new)));
	%max_vec_imag(ind) = max(max(imag(A_old - A_new)));
	%min_vec_real(ind) = min(min(real(A_old - A_new)));
	%min_vec_imag(ind) = min(min(imag(A_old - A_new)));
end
maximum
%figure
%hold all
%plot(wf, max_vec_real)
%plot(wf, min_vec_real)
%%xlim([-4 4])
%
%figure
%hold all
%plot(wf, max_vec_imag)
%plot(wf, min_vec_imag)
%%xlim([-4 4])

figure
hold all
plot(wf, diff)
%xlim([-100 100])

%Vergleiche feedback und central:

A_2mu = dy_new_data_short_str(3).m(pos_NfbP_2mu+1).m;

A_feedback = dy_new_data_long_str(3).m(1,1).m; 

diff_stat_p = max(max(abs(A_2mu - A_feedback)))



A_0 = dy_new_data_short_str(4).m(pos_NfbX_0+1).m;

A_feedback = dy_new_data_long_str(4).m(1,1).m; 

diff_stat_x = max(max(abs(A_0 - A_feedback)))

figure
imagesc(real(A_0 - A_feedback))

A_0 = dy_new_data_short_str(5).m(pos_NfbX_0+1).m;

A_feedback = dy_new_data_long_str(5).m(1,1).m; 

diff_stat_d = max(max(abs(A_0 - A_feedback)))
