PROGRAM=/home/hpc/uh3o1/ri26yad/bin/conductance_zero_mag_compare
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Conductance_zero_mag

#OPTIONS="/naslx/projects/uh3o1/ri26yad/DATA/keldysh_low_acc/Conductance/U_.55_B_0.00_T0/Chain_N31_Vg0.250000_mu-1.475000000_T0.000000000_h0.000000000_U\ 0.550000.mat"
OPTIONS="/naslx/projects/uh3o1/ri26yad/DATA/keldysh_low_acc/Conductance/U_.55_B_0.0_T_0.01_/Chain_N31_Vg0.250000_mu-1.475000000_T0.010000000_h0.000000000_U\ 0.550000.mat"
scf=bash_script
echo "#!/bin/bash" > $scf
cat SLURM_OPTIONS >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "$PROGRAM" "$OPTIONS" >> $scf
bash bash_script



