cmake_minimum_required(VERSION 3.10)

## Only one of these option can ever be set to ON! This is automatically done by the compile_cluster.sh script.
option(WORKSTATION "Build on workstation" ON)
option(ASC         "Build on ASC cluster" OFF)
option(KCS         "Build on KCS cluster" OFF)
option(SuperMUC    "Build on SuperMUC"    OFF)
option(JSC         "Build at JSC"         OFF)


## Set the C++ standard, the compiler, compiler optimization options (O0 - off, O3 - best) and add -fopenmp
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpicxx)

if(WORKSTATION)
    # for DEBUG build:
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -march=native -flto -fopenmp ")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O3 -march=native -flto ")
    # for RELEASE build:
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -flto -fopenmp")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -march=native -flto")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -march=native -flto")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -march=native -flto")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} -fopenmp")
endif()

if(WORKSTATION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} -fopenmp") # unnecessary??
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage") // for code coverage analysis (slows down the code a lot!)

    ## enable profiling (will generate a gmon.out file which can be interpreted with gprof in Linux).
    ## Comment out when no profiling shall be done!!
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

project(Keldysh_mfRG)

## include header and source code files
include_directories(${CMAKE_SOURCE_DIR})
file(GLOB_RECURSE keldysh_SRC *.c++)                                    # collects all filenames with extension .c++
file(GLOB_RECURSE test_SRC ${CMAKE_SOURCE_DIR}/tests/unit_tests/*.cxx)  # collects all filenames with extension .cxx
file(GLOB_RECURSE keldysh_H *.hpp *.h)                                  # collects all filenames with extension .h and .hpp


## standard executable for production runs
add_executable(Keldysh_mfRG main.cpp ${keldysh_H} ${keldysh_SRC})

## executable for (unit) tests (using Catch)
add_executable(Keldysh_mfRG_testing tests/unit_tests/unit_tests.cpp ${keldysh_H} ${keldysh_SRC} ${test_SRC})

## executable for production runs of FPP project
add_executable(FPP_mfRG FPP/FPP_main.cpp ${keldysh_H} ${keldysh_SRC})

## executable for integrand tests
add_executable(get_integrand_dGammaL tests/integrand_tests/get_integrand_dGammaL.cpp ${keldysh_H} ${keldysh_SRC})

## rest works so far only on workstation and ASC
if(WORKSTATION OR ASC)
    include_directories(/usr/include/c++/7/)
    include_directories(/usr/include/mpi/)
    include_directories(/usr/lib/gcc/x86_64-linux-gnu/7/include/) # omp.h
    include_directories(/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/include/)
    include_directories(/usr/include/eigen3 /usr/include/eigen3/unsupported)

    target_link_libraries(Keldysh_mfRG -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_mfRG -lgsl -lgslcblas)
    target_link_libraries(Keldysh_mfRG fftw3)

    target_link_libraries(Keldysh_mfRG_testing -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_mfRG_testing -lgsl -lgslcblas)
    target_link_libraries(Keldysh_mfRG_testing fftw3)

    target_link_libraries(FPP_mfRG -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(FPP_mfRG -lgsl -lgslcblas)
    target_link_libraries(FPP_mfRG fftw3)

    target_link_libraries(get_integrand_dGammaL -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(get_integrand_dGammaL -lgsl -lgslcblas)
    target_link_libraries(get_integrand_dGammaL fftw3)
elseif(KCS)
    include_directories($ENV{HDF5_INC} $ENV{GSL_INC} $ENV{BOOST_INC} $ENV{FFTW_INC} )
    include_directories($ENV{EIGEN_BASE}/include/eigen3)

    target_link_libraries(Keldysh_mfRG $ENV{HDF5_CPP_SHLIB} $ENV{HDF5_SHLIB} $ENV{SZIP_LIB} -lz)
    target_link_libraries(Keldysh_mfRG $ENV{FFTW_LIB} $ENV{GSL_LIB} -L$ENV{BOOST_LIBDIR}$)
elseif(JSC)
    include_directories(/p/software/juwelsbooster/stages/2022/software/HDF5/1.12.1-iompi-2021b/include) # $HDF5_INC
    include_directories(/p/software/juwelsbooster/stages/2022/software/FFTW/3.3.10-iompi-2021b/include) # $FFTW_INC
    include_directories(/p/software/juwelsbooster/stages/2022/software/GSL/2.7-GCCcore-11.2.0/include) # $GSL_INC
    include_directories(/p/software/juwelsbooster/stages/2022/software/Boost/1.78.0-GCCcore-11.2.0/include) # $BOOST_INC

    target_link_libraries(Keldysh_mfRG -L/p/software/juwelsbooster/stages/2022/software/HDF5/1.12.1-iompi-2021b/lib -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5) # HDF5
    target_link_libraries(Keldysh_mfRG -Wl,-static -L/p/software/juwelsbooster/stages/2022/software/FFTW/3.3.10-iompi-2021b/lib -lfftw3f -lfftw3 -lfftw3l -lm -Wl,-dy) # FFTW_LIB
    target_link_libraries(Keldysh_mfRG -L/p/software/juwelsbooster/stages/2022/software/GSL/2.7-GCCcore-11.2.0/lib -lgsl -lgslcblas) # GSL
    target_link_libraries(Keldysh_mfRG -L/p/software/juwelsbooster/stages/2022/software/Boost/1.78.0-GCCcore-11.2.0/lib$) # BOOST
endif()
