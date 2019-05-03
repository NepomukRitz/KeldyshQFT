PROGRAM=/home/hpc/uh3o1/ri26yad/bin/density_matrix
DIR_DATA_OUT=/naslx/projects/uh3o1/ri26yad/DATA/Density_matrix_eigenvalues/Without_leads_reg3

DIR_INPUT=/naslx/projects/uh3o1/ri26yad/DATA/Keldysh_long_range_equilibrium/Without_leads


for file in "$DIR_INPUT"/Chain_disorder_central_build_N61_Vg1.200000_Vsg0.000000_width0_mu0.000000000_T*_h0.000000000_U\ 0.800000_seed5.mat
do
    #echo "$file"
	OPTIONS=$(echo "$file" | sed -e 's/ /\\ /g')
	echo "$OPTIONS"
	scf=bash_script
	echo "#!/bin/bash" > $scf
	cat SLURM_OPTIONS >> $scf
	echo "cd $DIR_DATA_OUT" >> $scf 
	echo "$PROGRAM" "$OPTIONS" >> $scf
	bash bash_script

done


