PROGRAM=/home/hpc/uh3o1/ri26yad/bin/density_matrix_non_int
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/With_leads

N=30
Vg=1.2
mu=0.0
num_freq_pre=600000

for T in 0.100000000 0.158500000 0.251200000 0.398100000 0.631000000
do
 	echo $T
	OPTIONS="-N $N -Vg $Vg -mu $mu -T $T -num_freq_pre $num_freq_pre"
	scf=bash_script
	echo "#!/bin/bash" > $scf
	cat SLURM_OPTIONS >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "$PROGRAM $OPTIONS" >> $scf
	bash bash_script
done
