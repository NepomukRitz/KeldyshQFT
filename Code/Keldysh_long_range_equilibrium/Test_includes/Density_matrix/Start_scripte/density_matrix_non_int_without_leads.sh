PROGRAM=/home/hpc/uh3o1/ri26yad/bin/density_matrix_non_int
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/Without_leads
N=8
Vg=1.0
mu=0.0
num_freq_pre=600000
seed=1

for T in 0.1  
do
 	echo $T
	OPTIONS="-N $N -Vg $Vg -mu $mu -T $T -num_freq_pre $num_freq_pre -seed $seed"
	scf=bash_script
	echo "#!/bin/bash" > $scf
	cat SLURM_OPTIONS >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "$PROGRAM $OPTIONS" >> $scf
	bash bash_script
done
