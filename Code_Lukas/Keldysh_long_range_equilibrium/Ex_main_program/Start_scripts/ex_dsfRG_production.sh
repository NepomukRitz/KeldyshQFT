L=5
Lu=3
N=30
Nff=1500
NfbP=500
NfbX=500
num_freq_pre=30000
NL_full=5

Vg=0.25
h=0.0
mu_start=-1.50
mu_step=0.01
mu_end=-1.44
mulist=$(seq $mu_start $mu_step $mu_end)
T=0.0
U0=0.5
U1=0.3
Xi=5.0

accp=0.0001
accx=0.0001
accs=0.0001
tol=0.000001
Lambda_ini=100000.0
Lambda_fin=0.000000002

nodes=8


#runname="Schar_ex_dsfRG_L_${L}_Lu_${Lu}_N_${N}_Nff_${Nff}_NfbP_${NfbP}_NfbX_${NfbX}_num_freq_pre_${num_freq_pre}_NL_full_${NL_full}_Vg_${Vg}_h_${h}_mu_%f_T_${T}_U0_${U0}_U1_${U1}_Xi_${Xi}"

PROGRAM=/p/home/jusers/weidinger1/juwels/bin/ex_dsfRG_production_juwels
#DIR_DATA_OUT=/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp/
#DIR_LOG=/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp/AA_LOG
DIR_DATA_OUT=/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp_no_lr_extrapolation/
DIR_LOG=/p/scratch/chmu26/hmu261/Ex_DATA/QPC_zero_temp_no_lr_extrapolation/AA_LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

for mu in $mulist
do
	OPTIONS="-L $L -Lu $Lu -N $N -Nff $Nff -NfbP $NfbP -NfbX $NfbX -num_freq_pre $num_freq_pre -NL_full $NL_full -Vg $Vg -h $h -mu $mu -T $T -U0 $U0 -U1 $U1 -Xi $Xi -accp $accp -accx $accx -accs $accs -tol $tol -Lambda_ini $Lambda_ini -Lambda_fin $Lambda_fin" 
	scf=slurm_script
	echo "#!/bin/bash" > $scf
	echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
	echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
	echo "#SBATCH --nodes=$nodes" >> $scf
	echo "#SBATCH --ntasks=$nodes" >> $scf
	echo "#SBATCH --ntasks-per-node=1" >> $scf
	echo "#SBATCH --cpus-per-task=48" >> $scf
	echo "#SBATCH --time=05:59:00" >> $scf
	echo "#SBATCH --partition=batch" >> $scf
	echo "#SBATCH --account=hmu26" >> $scf
	echo "export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}" >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "srun $PROGRAM $OPTIONS" >> $scf
	sbatch slurm_script
done
