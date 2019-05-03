PROGRAM=/gpfs/homea/hmu26/hmu261/bin/conductance_provisorium
DIR_DATA_OUT=/gpfs/work/hmu26/hmu261/DATA/Unit_tests/Conductance_lead_system

OPTIONS="/gpfs/work/hmu26/hmu261/DATA/Final_test/dsfRG_mpi_L0_Lu0_N15_Nff_1500_NfbP_1500_NfbX_1500_num_freq_pre_30000_Vg_0.250000_h_0.000000_mu_-1.475000_T_0.002500_U0_0.500000_U1_0.000000_Xi_5.000000_Lambda_ini_99999.500001_Lambda_fin_0.000000_number_of_nodes_8.mat 0 0.5 0.0 5.0"
scf=bash_script
echo "#!/bin/bash" > $scf
#cat SLURM_OPTIONS >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "$PROGRAM" "$OPTIONS" >> $scf
bash bash_script



