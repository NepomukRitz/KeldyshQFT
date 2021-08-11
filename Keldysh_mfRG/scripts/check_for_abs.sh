#!/bin/bash

a=`grep -rn -o --include=*.h "abs(" ./ | wc -l`
grep -rn --include=*.h "abs(" ./ > a.txt
#echo $a

b=`grep -rn -o --include=*.h "std::abs(" ./ | wc -l`
grep -rn --include=*.h "std::abs(" ./ > b.txt  
#echo $b


c=`grep -rn -o --include=*.h "\.abs(" ./ | wc -l`
#echo $c

val=`expr $a - $b - $c - 2`





if [ $val == 0 ]
then
	echo "Passt."
	rm a.txt
	rm b.txt
else
	echo "Please check that all occuring 'abs' are from the vec class. (Otherwise you forgot a 'std::')"
	echo "There are " $val "too many 'abs()!'" 

	diff a.txt b.txt 
	
	rm a.txt
	rm b.txt

	exit 1
fi


