PROGRAM=/home/hpc/uh3o1/ri26yad/bin/density_matrix_non_int
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues
N=15
Vg=1.2
mu=0.0
T=0.0
num_freq_pre=600000
OPTIONS="-N $N -Vg $Vg -mu $mu -T $T -num_freq_pre $num_freq_pre"

scf=bash_script
echo "#!/bin/bash" > $scf
echo "cd $DIR_DATA_OUT" >> $scf 
echo "$PROGRAM $OPTIONS" >> $scf
bash bash_script
