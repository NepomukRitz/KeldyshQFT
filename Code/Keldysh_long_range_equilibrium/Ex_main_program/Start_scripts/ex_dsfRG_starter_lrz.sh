PROGRAM=~/bin/ex_dsfRG
DIR_DATA_OUT="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Ex_dsfRG/self_optimized"
DIR_LOG="/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/LOG"
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

scf=slurm_script
echo "#!/bin/bash" > $scf
echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
echo "#SBATCH -J MPI_job" >> $scf
echo "#SBATCH --get-user-env" >> $scf
echo "#SBATCH --clusters=mpp2" >> $scf
echo "#SBATCH --nodes=1" >> $scf
echo "#SBATCH --ntasks=1" >> $scf
echo "#SBATCH --ntasks-per-node=1" >> $scf
echo "#SBATCH --cpus-per-task=28" >> $scf
echo "#SBATCH --time=00:15:00" >> $scf
echo "export OMP_NUM_THREADS=\${SLURM_CPUS_PER_TASK}" >> $scf
echo "source /etc/profile.d/modules.sh" >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "mpiexec $PROGRAM" >> $scf
sbatch slurm_script
