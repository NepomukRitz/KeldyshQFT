NODES=5
L=1
Lu=0
N=4
Nff=1500
NfbP=1500
NfbX=1500
num_freq_pre=30000
Vg=0.25
h=0.0
T=0.0025
U0=0.8
U1=0.0
Xi=5.0
mu_start=-1.475
mu_step=0.01
mu_end=-1.475
mulist=$(seq $mu_start $mu_step $mu_end)

runname="Schar_dsfRG_mpi_L${L}_Lu${Lu}_N${N}_Nff${Nff}_NfbP${NfbP}_NfbX${NfbX}_num_freq_pre${num_freq_pre}_Vg${Vg}_h${h}_T${T}_U0${U0}_U1${U1}_Xi${Xi}_mu${mu_start}_${mu_step}_${mu_end}_nodes${NODES}"
PROGRAM=/gpfs/homea/hmu26/hmu261/bin/dsfRG_mpi
DIR_DATA_OUT=/gpfs/work/hmu26/hmu261/DATA/Static_long_range/QPC/Test_complete_MPI/$runname
DIR_LOG=/gpfs/work/hmu26/hmu261/DATA/Static_long_range/QPC/Test_complete_MPI/$runname/LOG
#DIR_DATA_OUT=/work/hmu26/hmu261/DATA/Static_long_range/QPC/Long_range_convergence/$runname
#DIR_LOG=/work/hmu26/hmu261/DATA/Static_long_range/QPC/Long_range_convergence/$runname/LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

for mu in $mulist
do
	OPTIONS="-L $L -Lu $Lu -N $N -Nff $Nff -NfbP $NfbP -NfbX $NfbX -num_freq_pre $num_freq_pre -Vg $Vg -h $h -mu $mu -T $T -U0 $U0 -U1 $U1 -Xi $Xi" 
	scf=slurm_script
	echo "#!/bin/bash" > $scf
	echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
	echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
	echo "#SBATCH --nodes=$NODES" >> $scf
	echo "#SBATCH --ntasks=$NODES" >> $scf
	echo "#SBATCH --cpus-per-task=24" >> $scf
	echo "#SBATCH --time=00:15:00" >> $scf
	echo "#SBATCH --partition=batch" >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "srun $PROGRAM $OPTIONS" >> $scf
	sbatch slurm_script
done