clear all
close all
load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_evaluate_right_hand_side_rpa_L_5_Lu_5_N_10_Nff_200_NfbP_200_NfbX_200_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.400000_U1_0.400000_Xi_5.000000_T_0.030000_Lambda_0.010000000.mat')

for ind=1:length(wbP)
 	A(ind) = dgamma_data_short_str(3).m(ind).m(11,11)/Measure_Flow;
end

figure
hold all
plot(wbP,real(A))
plot(wbP,imag(A))
xlim([-4 4])

%load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_intertwined_L_5_Lu_5_N_10_Nff_200_NfbP_200_NfbX_200_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.400000_U1_0.400000_Xi_5.000000_T_0.030000_Lambda_0.010000.mat')
load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_L_5_Lu_5_N_10_Nff_200_NfbP_200_NfbX_200_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.400000_U1_0.400000_Xi_5.000000_T_0.030000_Lambda_0.009999900.mat')

for ind=1:length(wbP)
 	B(ind) = gamma_data_short_str(3).m(ind).m(11,11);
end



%load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_intertwined_L_5_Lu_5_N_10_Nff_200_NfbP_200_NfbX_200_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.400000_U1_0.400000_Xi_5.000000_T_0.030000_Lambda_0.010010.mat')
load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/P_rpa_zero_mag/P_rpa_zero_mag_L_5_Lu_5_N_10_Nff_200_NfbP_200_NfbX_200_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_U0_0.400000_U1_0.400000_Xi_5.000000_T_0.030000_Lambda_0.010000100.mat')

for ind=1:length(wbP)
 	C(ind) = gamma_data_short_str(3).m(ind).m(11,11);
end
dA = (C-B)/2e-7;
plot(wbP, real(dA))
plot(wbP, imag(dA))
Measure_Flow

[diff a] = max(abs(A -dA))
A(a)
dA(a)
wbP(a)
