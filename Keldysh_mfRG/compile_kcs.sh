#!/bin/bash

module load hdf5/serial
module load fftw

export LANG=C
export LC_ALL=C

HDF5="$HDF5_INC $HDF5_CPP_LIB $HDF5_LIB $SZIP_LIB -lz"
FFTW="$FFTW_INC $FFTW_LIB"

mpiCC --std=c++11  main.cpp -o main.o -fopenmp $HDF5 $FFTW
