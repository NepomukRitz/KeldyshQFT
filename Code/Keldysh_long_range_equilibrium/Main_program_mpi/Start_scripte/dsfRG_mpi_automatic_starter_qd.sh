NODES=5
L=20
Lu=0
N=40
Nff=1500
NfbP=1500
NfbX=1500
num_freq_pre=30000
mu=-1.0
V_sg=1.007
pot_width=20
h=0.0
T=0.005
U0=1.0
U1=0.0
Xi=5.0
Vg_start=0.6
Vg_step=0.01
Vg_end=1.0
Vglist=$(seq $Vg_start $Vg_step $Vg_end)

runname="Schar_dsfRG_mpi_L${L}_Lu${Lu}_N${N}_Nff${Nff}_NfbP${NfbP}_NfbX${NfbX}_num_freq_pre${num_freq_pre}_mu${mu}_h${h}_T${T}_U0${U0}_U1${U1}_Xi${Xi}_Vsg${Vsg}_pot_width${pot_width}_Vg${Vg_start}_${Vg_step}_${Vg_end}_nodes${NODES}"
PROGRAM=/homea/hmu26/hmu261/bin/dsfRG_mpi
#DIR_DATA_OUT=/work/hmu26/hmu261/DATA/Static_long_range/QPC/Test_complete_MPI/$runname
#DIR_LOG=/work/hmu26/hmu261/DATA/Static_long_range/QPC/Test_complete_MPI/$runname/LOG
DIR_DATA_OUT=/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/$runname
DIR_LOG=/work/hmu26/hmu261/DATA/Static_long_range/QD/Long_range_convergence/$runname/LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

for Vg in $Vglist
do
	OPTIONS="-L $L -Lu $Lu -N $N -Nff $Nff -NfbP $NfbP -NfbX $NfbX -num_freq_pre $num_freq_pre -Vg $Vg -h $h -mu $mu -T $T -U0 $U0 -U1 $U1 -Xi $Xi -V_sg $V_sg -pot_width $pot_width" 
	scf=slurm_script
	echo "#!/bin/bash" > $scf
	echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
	echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
	echo "#SBATCH --nodes=$NODES" >> $scf
	echo "#SBATCH --ntasks=$NODES" >> $scf
	echo "#SBATCH --cpus-per-task=24" >> $scf
	echo "#SBATCH --time=08:00:00" >> $scf
	echo "#SBATCH --partition=batch" >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "srun $PROGRAM $OPTIONS" >> $scf
	sbatch slurm_script
done
