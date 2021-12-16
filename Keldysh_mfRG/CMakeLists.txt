cmake_minimum_required(VERSION 3.10)
project(Keldysh_mfRG)

set(CMAKE_CXX_STANDARD 17)

## Set optimization level of compiler:
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

set(CMAKE_CXX_COMPILER mpicxx)

## added -fopenmp
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} -fopenmp") # unnecessary??
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage") // for code coverage analysis (slows down the code a lot!)

#set(SOURCE_FILES main.cpp) # diagrammatic_combinations.h interpolations.h testFunctions.h correctionFunctions.h)
#add_executable(openmp_test_clion ${SOURCE_FILES})

file(GLOB keldysh_SRC CONFIGURE_DEPENDS "*.h" "*.hpp" "*.c++")  # collects all filenames with extension .h, .hpp and .c++

## standard executable for integration tests and production runs
add_executable(Keldysh_mfRG
        main.cpp
        ${keldysh_SRC} grids/momentum_grid.c++ symmetries/symmetry_transformations.c++ postprocessing/KramersKronig.c++ utilities/write_data2file.c++ mfRG_flow/flow.c++ utilities/hdf5_routines.c++)

## executable for unit tests using Catch
add_executable(Keldysh_mfRG_testing
        tests/unit_tests/unit_tests.cpp
        ${keldysh_SRC} grids/momentum_grid.c++ symmetries/symmetry_transformations.c++ postprocessing/KramersKronig.c++ utilities/write_data2file.c++ mfRG_flow/flow.c++ utilities/hdf5_routines.c++)

## standard executable for integration tests and production runs
add_executable(FPP_mfRG
        FPP/FPP_main.cpp
        ${keldysh_SRC} grids/momentum_grid.c++ symmetries/symmetry_transformations.c++ postprocessing/KramersKronig.c++ postprocessing/postprocessing.c++ utilities/write_data2file.c++ mfRG_flow/flow.c++ utilities/hdf5_routines.c++ tests/test_PrecalculatedBubble.c++ tests/test_symmetries.c++)

add_executable(get_integrand_dGammaL
        tests/integrand_tests/get_integrand_dGammaL.cpp
        ${keldysh_SRC} grids/momentum_grid.c++ symmetries/symmetry_transformations.c++ postprocessing/KramersKronig.c++ utilities/write_data2file.c++ mfRG_flow/flow.c++ utilities/hdf5_routines.c++ tests/integrand_tests/saveIntegrand.c++)

include_directories(/usr/include/mpi/)
include_directories(/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/include)

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