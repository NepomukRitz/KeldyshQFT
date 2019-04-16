DIR_CURRENT=$PWD

foldername=$1
cd $foldername

for filename in X*.mat; do
	echo $filename
	condname="Conductance/cond_$filename"
	if [ -d $condname ]
	then
		echo "$condname already exists"
	else
        	/p/home/jusers/weidinger1/juwels/bin/Compute_Conductance $foldername $filename 
	fi
done

mv cond* ./Conductance

cd $DIR_CURRENT

