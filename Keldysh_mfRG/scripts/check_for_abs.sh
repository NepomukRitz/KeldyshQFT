#!/bin/bash
### This script makes sure that only 'std::abs()' is used for scalars throughout the code ('std::abs' is GOOD abs).
### 'std::fabs()' and 'abs()' are to be avoided (these are BAD abs). (Otherwise the compiler runs into problems with choosing the correct definition for this function.)
### Only for vector is a allowed to use something else: e.g. foo.abs()       (These are also GOOD abs)

## Search for occurences of 'abs(' (This finds both GOOD abs and BAD abs, i.e. 'std::abs(', 'fabs(', '.abs(')
a=`grep -rn -o --include=*.h "abs(" ./ | wc -l`
## write the occurrences in a file
grep -rn --include=*.h "abs(" ./ > a.txt

## Search for occurences of 'std::abs(' (GOOD abs)
b=`grep -rn -o --include=*.h "std::abs(" ./ | wc -l`
## write the occurrences in a file
grep -rn --include=*.h "std::abs(" ./ > b.txt  

## Search for occurences of '.abs(' (GOOD abs)
c=`grep -rn -o --include=*.h "\.abs(" ./ | wc -l`

## Compute the number of BAD abs; the 2 come from the declaration and definition of vec<double>::abs() and are hence GOOD abs
val=`expr $a - $b - $c - 2`





if [ $val == 0 ]
then
	echo "Good! We only use 'std::abs()' for scalars."
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


