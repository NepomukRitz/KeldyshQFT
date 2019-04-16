PROGRAM=/gpfs/homea/hmu26/hmu261/bin/dsfRG_mpi
DIR_DATA_OUT=/gpfs/work/hmu26/hmu261/DATA/Final_test
DIR_LOG=/gpfs/work/hmu26/hmu261/DATA/Test_LOG
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

scf=slurm_script
echo "#!/bin/bash" > $scf
echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
cat SLURM_OPTIONS_REST >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "srun $PROGRAM" >> $scf
sbatch slurm_script
