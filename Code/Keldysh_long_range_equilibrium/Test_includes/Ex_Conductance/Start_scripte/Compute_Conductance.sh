DIR_CURRENT=$PWD
DIR_OUTPUT="/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_Conductance"
#DIR_OUTPUT="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Conductance"
cd $DIR_OUTPUT
#~/bin/Compute_Conductance "/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_Conductance" "X_L0_Lu0_N30_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu-1.5000_T0.005000_Uc0.65_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
~/bin/Compute_Conductance "/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation_old2" "X_L0_Lu0_N15_Nff1500_NfbP1500_NfbX1500_pre30000_NL0_Vg0.2500_h0.000000_mu-1.4750_T0.010000_Uc0.55_Uo0.00_Xi5.00_ap1e-04_ax1e-04_as1e-04_tol1e-06_Li100000.0_Lf0.0_no1.mat"
cd $DIR_CURRENT

