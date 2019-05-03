PROGRAM=/gpfs/homea/hmu26/hmu261/bin/conductance_provisorium
DIR_DATA_IN=/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0/Schar_dsfRG_mpi_L10_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/
DIR_DATA_OUT=/gpfs/work/hmu26/hmu261/DATA_PRODUCTION/QPC_short_interactions/T0.0_conductance/Conductance_schar_dsfRG_mpi_L10_Lu0_N30_Nff1500_NfbP1500_NfbX1500_num_freq_pre30000_Vg0.25_h0.0_T0.0_U00.65_U10.0_Xi5.0_mu-1.6_0.005_-1.4_nodes8/

Lu=0
U0=0.65
U1=0.0
Xi=5.0
mkdir -p $DIR_DATA_OUT

for filename in $DIR_DATA_IN*.mat; do
	OPTIONS="$filename $Lu $U0 $U1 $Xi"
	scf=bash_script
	echo "#!/bin/bash" > $scf
	#cat SLURM_OPTIONS >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "$PROGRAM" "$OPTIONS" >> $scf
	bash bash_script
done




