#!/bin/bash

# Shell script for the compilation process on the ASC-cluster, the KCS-cluster and the SuperMUC
# Usage: ...

# check option and load necessary modules
if [ "$1" == "--ASC" ]
  then
  export CLUSTER=ASC
  module load cmake
elif [ "$1" == "--KCS" ]
then
  export CLUSTER=KCS
  module unload spack gcc hdf5 fftw gsl boost cmake # in case some old versions have been loaded previously
  module unload devEnv/Intel/2019 itac/2019

  module load spack/22.2.1
  module load intel-mpi/2019-intel
  module load gcc/9
  #module load hdf5/1.10.7-intel21-impi
  module load hdf5/1.8.22-gcc11
  module load fftw/3.3.10
  module load gsl/2.7-intel21
  module load boost/1.75.0-intel21-impi
  module load eigen/3.4.0-intel21
  module load cmake/3.21.4

  #module load intel intel-mpi/2019.8.254

  #module load gcc/9
  #module load hdf5/1.8.21-gcc8
  #module load fftw
  #module load gsl
  #module load boost/1.70.0-intel19
  #module load cmake
  #module load eigen/3.3.7-intel19
elif [ "$1" == "--JSC" ]
then
  export CLUSTER=JSC
  module load Stages/2022 Intel/2021.4.0 OpenMPI/4.1.2
  module load HDF5/1.12.1 # could be problematic, because already version 1.10 (instead of 1.08) and no threadsafe option
  module load FFTW
  module load GSL
  module load Boost
  module load CMake/3.21.1
else
  echo "Invalid cluster option! Needs to be --ASC, --KCS or --JSC"
  exit
fi

export LANG=C
export LC_ALL=C

dir_path=$(dirname "$(realpath "$0")")    # directory of this script
src_dir="$dir_path"/..                    # source directory
build_dir="$src_dir"/build                # build directory
mkdir "$build_dir"

# navigate to the directory with the CMakeLists.txt file
cd "$build_dir" || exit


# load cmake project
cmake -DWORKSTATION=OFF -D$CLUSTER=ON -S$src_dir -B$build_dir


# build
TARGET=source
cmake --build . --target $TARGET -- -j 9

mv ./$TARGET $src_dir

GREEN='\033[1;32m' # green
NC='\033[0m' # no color
echo -e "\n Moved executable ${GREEN}$TARGET ${NC}to Keldysh_mfRG source directory.\n"
