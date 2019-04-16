PROGRAM=/p/home/jusers/weidinger1/juwels/bin/ex_dsfRG_vgl_flow_old
DIR_DATA_OUT=/p/scratch/chmu26/hmu261/Ex_Unit_tests/Ex_dsfRG
DIR_LOG=/p/scratch/chmu26/hmu261/Ex_Unit_tests/LOG
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
