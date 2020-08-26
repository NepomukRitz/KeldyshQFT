#!/bin/bash

module load hdf5/1.10.5-gcc

export LANG=C
export LC_ALL=C

MPI="-I/usr/include/mpi/"
HDF5="-I/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/include -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp"
FFTW="-lfftw3"

mpiCC --std=c++11  main.cpp -o main.o -fopenmp $MPI $HDF5 $FFTW

