clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_vgl_flow_old.mat')

%Vgl short data:

for ind=1:length(wf)
 	Self_up_diff(ind) = max(max(abs( y_old(1).m(ind).m - gamma_data_short_str(8).m(ind).m )));
 	Self_down_diff(ind) = max(max(abs( y_old(2).m(ind).m - gamma_data_short_str(9).m(ind).m )));
 	Pud_diff(ind) = max(max(abs( y_old(3).m(ind).m - U/4 - gamma_data_short_str(3).m(ind).m )));
 	Xud_diff(ind) = max(max(abs( y_old(4).m(ind).m - U/4 - gamma_data_short_str(4).m(ind).m )));
 	Duu_diff(ind) = max(max(abs( y_old(5).m(ind).m - gamma_data_short_str(5).m(ind).m )));
 	Ddd_diff(ind) = max(max(abs( y_old(6).m(ind).m - gamma_data_short_str(6).m(ind).m )));
end


for ind1=1:length(gamma_data_long_str(3).m)
	for ind2=1:length(gamma_data_long_str(3).m)
 		Puu_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(1).m(ind1,ind2).m ) ) ); 
 		Pdd_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(2).m(ind1,ind2).m ) ) );
 		Dud_long_max(ind1,ind2) = max(max( abs( gamma_data_long_str(7).m(ind1,ind2).m ) ) );
	end
end

