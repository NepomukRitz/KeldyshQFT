DIR_CURRENT=$PWD

foldername=$1
cd $foldername

for filename in X*.mat; do
	echo $filename
	condname="Conductance_old/cond_$filename"
	if [ -d $condname ]
	then
		echo "$condname already exists"
	else
			/home/hpc/uh3o1/ri26yad/bin/Compute_Conductance $foldername $filename
			mv cond* ./Conductance_old
	fi
done


cd $DIR_CURRENT

