clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG.mat')

for ind=1:length(gamma_data_short_str(8).m)
 	B(ind) = gamma_data_short_str(8).m(ind).m(21,21);
 	C(ind) = gamma_data_before_flow_short_str(8).m(ind).m(21,21);
end

max(abs(B - C))

figure
hold 
plot(wf,real(B-C))
plot(wf,imag(B-C))

