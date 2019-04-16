clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_vgl_flow_old.mat')

for ind=1:length(y_old(3).m)
 	A(ind) = y_old(1).m(ind).m(11,11);
 	%A(ind) = y_old(3).m(ind).m(11,11) - U(11,11)/4;
end

for ind=1:length(gamma_data_short_str(3).m)
 	B(ind) = gamma_data_short_str(8).m(ind).m(11,11);
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

diff = 0.0;
for ind=1:length(y_old(1).m)
 	A = y_old(1).m(ind).m; 
 	B = gamma_data_short_str(8).m(ind).m;
	diff = max(diff, max(max(abs(A-B))));
end
diff

