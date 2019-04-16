clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_flow_new.mat')


for ind=1:length(gamma_data_short_str(3).m)
 	B(ind) = gamma_data_short_str(3).m(ind).m(11,11);
end


figure
plot(wf,imag(B))
%xlim([-10 10])

figure
plot(wf,real(B))
%xlim([-10 10])



