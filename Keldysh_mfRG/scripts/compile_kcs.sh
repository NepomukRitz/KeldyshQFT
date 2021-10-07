#!/bin/bash
#  environment variable KELDYSH_MFRG needs to point to the "Keldysh_mfRG" directory of the repository:
#  in ~/.bashrc:
#  export KELDYSH_MFRG="/dss/dsskcsfs01/pn34vu/pn34vu-dss-0000/<lrz-ID>/mfrg/Keldysh_mfRG"

module load gcc
module load hdf5/1.8.20-cxx-frt-threadsafe
module load fftw
module load gsl
module load boost/1.61_icc

export LANG=C
export LC_ALL=C

HDF5="$HDF5_INC $HDF5_CPP_SHLIB $HDF5_SHLIB $SZIP_LIB -lz"
FFTW="$FFTW_INC $FFTW_LIB"
GSL="$GSL_INC $GSL_LIB"
BOOST="$BOOST_INC -L$BOOST_LIBDIR$"

mpiCC -std=c++17 $KELDYSH_MFRG/main.cpp -o $KELDYSH_MFRG/main.o -fopenmp $FFTW $HDF5 $GSL $BOOST


