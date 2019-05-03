clear all
close all

load('/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_simple/dsfRG_vgl_right_hand_side.mat')

%X_dyn = dy_old(4).m(pos_NfbX_0+1).m;
X_dyn = dy_old(4).m(752).m;
X_stat = dy_new_data_long_str(4).m(1,1).m;

x= max(max(abs(X_dyn - X_stat )))


%D_dyn = dy_old(5).m(pos_NfbX_0+1).m;
D_dyn = dy_old(5).m(752).m;
D_stat = dy_new_data_long_str(5).m(1,1).m;

d= max(max(abs(D_dyn - D_stat )))

P_dyn = dy_old(3).m(382).m;
P_stat = dy_new_data_long_str(3).m(1,1).m;

p= max(max(abs(P_dyn - P_stat )))

figure
imagesc(real(P_dyn - P_stat))
figure
imagesc(imag(P_dyn - P_stat))
