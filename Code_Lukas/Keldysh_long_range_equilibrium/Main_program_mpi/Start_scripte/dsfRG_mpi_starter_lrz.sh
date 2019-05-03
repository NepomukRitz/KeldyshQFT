PROGRAM=/home/hpc/uh3o1/ri26yad/bin/dsfRG_mpi
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/Main_program_mpi
DIR_LOG=/home/hpc/uh3o1/ri26yad/LOG/Main_program_mpi
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 

scf=slurm_script
echo "#!/bin/bash" > $scf
echo "#SBATCH --output=$DIR_LOG/mpi-out.%j" >> $scf
echo "#SBATCH --error=$DIR_LOG/mpi-err.%j" >> $scf
cat SLURM_OPTIONS_LRZ_REST >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "srun $PROGRAM" >> $scf
sbatch slurm_script
