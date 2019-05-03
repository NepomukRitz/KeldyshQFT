PROGRAM=/home/hpc/uh3o1/ri26yad/bin/hf_equations
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Unit_tests/HF_equations
DIR_LOG=/home/hpc/uh3o1/ri26yad/LOG/LOG_HF_equations
mkdir -p $DIR_LOG 
mkdir -p $DIR_DATA_OUT 
scf=slurm_script
echo "#!/bin/bash" > $scf
echo "#SBATCH -o $DIR_LOG/myjob.%j.%N.out" >> $scf
cat SLURM_OPTIONS_REST >> $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "$PROGRAM" >> $scf
sbatch slurm_script
