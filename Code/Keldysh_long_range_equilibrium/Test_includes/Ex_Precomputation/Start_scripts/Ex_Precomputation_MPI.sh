PROGRAM=/p/home/jusers/weidinger1/juwels/bin/Ex_Precomputation
DIR_DATA_OUT=/p/scratch/chmu26/hmu261/DATA/Unit_tests/Ex_Precomputation
DIR_LOG=/p/scratch/chmu26/hmu261/DATA/Unit_tests/Ex_Precomputation/LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

scf=slurm_script
echo "#!/bin/bash" > $scf
echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
cat SLURM_OPTIONS_REST >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "srun $PROGRAM" >> $scf
sbatch -A hmu26 slurm_script

