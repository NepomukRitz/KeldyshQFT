clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/check_precomputation_in_flow.mat')


for ind=1:length(external_frequencies_Gu_new)
 	A_new(ind) = external_frequencies_Gu_new(ind).m(31,31);
 	A_old(ind) = external_frequencies_Gu_old(ind).m(31,31);
end

figure
hold 
plot(external_frequencies_new, imag(A_new))
plot(external_frequencies_old, imag(A_old))
xlim([-4 4])

max(abs(A_new -A_old))
