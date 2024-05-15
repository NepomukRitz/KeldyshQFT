cmake_minimum_required(VERSION 3.10)

## Only one of these option can ever be set to ON! This is automatically done by the compile_cluster.sh script.
option(WORKSTATION "Build on workstation" OFF)
option(ASC         "Build on ASC cluster" OFF)
option(KCS         "Build on KCS cluster" OFF)
option(SuperMUC    "Build on SuperMUC"    OFF)
option(JSC         "Build at JSC"         OFF)

option(MACBOOK     "Build on MacBook"     ON)


## Set the C++ standard, the compiler, compiler optimization options (O0 - off, O3 - best) and add -fopenmp
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(WORKSTATION)
    # for DEBUG build:
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O3 -march=native -flto  -fopenmp ")
    set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -D_DEBUG -O3 -march=native -flto  ")
    # for RELEASE build:
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ")
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
file(GLOB_RECURSE test_SRC ${CMAKE_SOURCE_DIR}/tests/unit_tests/*.cxx ${CMAKE_SOURCE_DIR}/tests/*.cxx)  # collects all filenames with extension .cxx
file(GLOB_RECURSE keldysh_H *.hpp *.h)                                  # collects all filenames with extension .h and .hpp


## standard executable for production runs
add_executable(Keldysh_mfRG main.cpp ${keldysh_H} ${keldysh_SRC})

## executable for post-processing
add_executable(Keldysh_postproc postprocessing/main_postproc.cpp ${keldysh_H} ${keldysh_SRC})

## executable for (unit) tests (using Catch)
add_executable(Keldysh_mfRG_testing tests/unit_tests/unit_tests.cpp ${keldysh_H} ${keldysh_SRC} ${test_SRC})


## executable for integrand tests
add_executable(get_integrand_dGammaL tests/integrand_tests/get_integrand_dGammaL.cpp ${keldysh_H} ${keldysh_SRC})
add_executable(get_integrand_dGamma_1loop tests/integrand_tests/get_integrand_dGamma_1loop.cpp ${keldysh_H} ${keldysh_SRC})
add_executable(get_integrand_Gamma1_Pi_Gamma2 tests/integrand_tests/get_integrand_Gamma1_Pi_Gamma2.cpp ${keldysh_H} ${keldysh_SRC})
add_executable(get_integrand_Sigma_SDE tests/integrand_tests/get_integrand_Sigma_SDE.cpp ${keldysh_H} ${keldysh_SRC})

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

    target_link_libraries(Keldysh_postproc -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_postproc -lgsl -lgslcblas)
    target_link_libraries(Keldysh_postproc fftw3)

    target_link_libraries(Keldysh_mfRG_testing -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_mfRG_testing -lgsl -lgslcblas)
    target_link_libraries(Keldysh_mfRG_testing fftw3)

    target_link_libraries(get_integrand_dGammaL -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(get_integrand_dGammaL -lgsl -lgslcblas)
    target_link_libraries(get_integrand_dGammaL fftw3)
    target_link_libraries(get_integrand_Gamma1_Pi_Gamma2 -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(get_integrand_Gamma1_Pi_Gamma2 -lgsl -lgslcblas)
    target_link_libraries(get_integrand_Gamma1_Pi_Gamma2 fftw3)
    target_link_libraries(get_integrand_dGamma_1loop -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(get_integrand_dGamma_1loop -lgsl -lgslcblas)
    target_link_libraries(get_integrand_dGamma_1loop fftw3)
    target_link_libraries(get_integrand_Sigma_SDE -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(get_integrand_Sigma_SDE -lgsl -lgslcblas)
    target_link_libraries(get_integrand_Sigma_SDE fftw3)
elseif(KCS)
    include_directories($ENV{HDF5_BASE}/include $ENV{GSL_BASE}/include $ENV{BOOST_INCDIR} $ENV{FFTW_BASE}/include)
    include_directories($ENV{EIGEN_BASE}/include/eigen3)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} $ENV{HDF5_CPP_SHLIB} $ENV{HDF5_SHLIB}") # sets necessary compiler flags for HDF5 library

    target_link_libraries(Keldysh_mfRG $ENV{SZIP_LIB} -lz)
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
elseif(MACBOOK)
    include_directories(
            /usr/local/include/
            /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/
            /opt/homebrew/Cellar/gsl/2.7.1/include
            /opt/homebrew/Cellar/hdf5/1.14.3/include
            /opt/homebrew/Cellar/boost/1.84.0/include
    )
    target_link_libraries(Keldysh_mfRG -L/opt/homebrew/Cellar/hdf5/1.14.3/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_mfRG -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas)

    target_link_libraries(Keldysh_postproc -L/opt/homebrew/Cellar/hdf5/1.14.3/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_postproc -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas)

    target_link_libraries(Keldysh_mfRG_testing -L/opt/homebrew/Cellar/hdf5/1.14.3/lib -lhdf5 -lhdf5_cpp)
    target_link_libraries(Keldysh_mfRG_testing -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas)
endif()


# Documentation
find_package(Doxygen REQUIRED)
get_target_property(KELDYSH_MFRG_PUBLIC_HEADER_DIR Keldysh_mfRG INTERFACE_INCLUDE_DIRECTORIES)
file(GLOB_RECURSE KELDYSH_MFRG_PUBLIC_HEADERS ${KELDYSH_MFRG_PUBLIC_HEADER_DIR}/*.h)
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/../Documentation)
#This will be the main output of our command
set(DOXYGEN_INDEX_FILE ${CMAKE_CURRENT_BINARY_DIR}/../Documentation/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${KELDYSH_MFRG_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/Documentation)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/../../docs)

add_custom_target(Sphinx ALL
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
                  # Tell Breathe where to find the Doxygen output
                  -Dbreathe_projects.KeldyshDiag=${DOXYGEN_OUTPUT_DIR}/xml/
                  ${SPHINX_SOURCE} ${SPHINX_BUILD}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  COMMENT "Generating documentation with Sphinx")