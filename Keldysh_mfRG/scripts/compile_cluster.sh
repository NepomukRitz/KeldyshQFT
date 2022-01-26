#!/bin/bash

# Shell script for the compilation process on the ASC-cluster, the KCS-cluster and the SuperMUC
# Usage: ...

# check option and load necessary modules
if [ "$1" == "--ASC" ]
  then
  export CLUSTER=ASC
elif [ "$1" == "--KCS" ]
then
  export CLUSTER=KCS
  module load gcc
  module load hdf5/1.8.20-cxx-frt-threadsafe
  module load fftw
  module load gsl
  module load boost/1.61_icc
else
  echo "Invalid cluster option! Needs to be --ASC or --KCS."
  exit
fi
module load cmake


export LANG=C
export LC_ALL=C


# navigate to the directory with the CMakeLists.txt file
cd ../


# load cmake project
cmake -DWORKSTATION=OFF -D$CLUSTER=ON .


# build
cmake --build . --target Keldysh_mfRG