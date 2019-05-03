#!/bin/bash

PATH_TO_LIB=$(pwd)/Code/lib_without_matlab
PATH_TO_EIGEN=$(pwd)/eigen

FILENAME=SIAM_constant_vertex

# compile using Lukas' "matrix" class
# g++ -std=c++11 -O0 -o $FILENAME $FILENAME.cpp $PATH_TO_LIB/matrix.cpp $PATH_TO_LIB/blnla.cpp  -I $PATH_TO_LIB

# compile using "eigen" library // TODO: "eigen" probably unnecessary
g++ -std=c++11 -O3 -o $FILENAME $FILENAME.cpp -I $PATH_TO_EIGEN

# execute file
./$FILENAME
