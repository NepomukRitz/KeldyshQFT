L=5
Lu=3
N=30
Nff=1500
NfbP=1500
NfbX=1500
num_freq_pre=30000
NL_full=0

Vg=0.25
h=0.0
mu_start=-1.55
mu_step=0.01
mu_end=-1.40
mulist=$(seq $mu_start $mu_step $mu_end)
T=0.0025
U0=0.5
U1=0.3
Xi=5.0

accp=0.0001
accx=0.0001
accs=0.0001
tol=0.000001
Lambda_ini=100000.0
Lambda_fin=0.000000002

#runname="Schar_ex_dsfRG_L_${L}_Lu_${Lu}_N_${N}_Nff_${Nff}_NfbP_${NfbP}_NfbX_${NfbX}_num_freq_pre_${num_freq_pre}_NL_full_${NL_full}_Vg_${Vg}_h_${h}_mu_%f_T_${T}_U0_${U0}_U1_${U1}_Xi_${Xi}"

PROGRAM=/home/hpc/uh3o1/ri26yad/bin/ex_dsfRG_production
#DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/Ex_DATA_lrz/QPC_zero_tmp_lrz/
#DIR_LOG=/naslx/projects/uh3o1/ri26yad/Ex_DATA_lrz/QPC_zero_tmp_lrz/AA_LOG
#DIR_DATA_OUT=/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/
#DIR_LOG=/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_zero_tmp_lrz_no_lr_extrapolation/AA_LOG
DIR_DATA_OUT=/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation_with_additional_freq_in_px/
DIR_LOG=/gpfs/scratch/uh3o1/ri26yad/Ex_DATA_lrz_scratch/QPC_finite_tmp_lrz_no_lr_extrapolation_with_additional_freq_in_px/AA_LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

for mu in $mulist
do
	OPTIONS="-L $L -Lu $Lu -N $N -Nff $Nff -NfbP $NfbP -NfbX $NfbX -num_freq_pre $num_freq_pre -NL_full $NL_full -Vg $Vg -h $h -mu $mu -T $T -U0 $U0 -U1 $U1 -Xi $Xi -accp $accp -accx $accx -accs $accs -tol $tol -Lambda_ini $Lambda_ini -Lambda_fin $Lambda_fin" 
	scf=slurm_script
	echo "#!/bin/bash" > $scf
	echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
	echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
	echo "#SBATCH -J MPI_job" >> $scf
	echo "#SBATCH --get-user-env" >> $scf
	echo "#SBATCH --clusters=mpp2" >> $scf
	echo "#SBATCH --nodes=4" >> $scf
	echo "#SBATCH --ntasks=4" >> $scf
	echo "#SBATCH --ntasks-per-node=1" >> $scf
	echo "#SBATCH --cpus-per-task=28" >> $scf
	echo "#SBATCH --time=05:59:59" >> $scf
	echo "export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}" >> $scf
	echo "source /etc/profile.d/modules.sh" >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "mpiexec $PROGRAM $OPTIONS" >> $scf
	sbatch slurm_script
done
