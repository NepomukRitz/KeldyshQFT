cmake_minimum_required(VERSION 3.14)
project(Keldysh_mfRG)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_COMPILER mpicxx)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} -std=c++14 -fopenmp")

#include_directories("$ gsl-2.5/")

set(SOURCE_FILES main.cpp diagrammatic_combinations.h interpolations.h testFunctions.h)
add_executable(openmp_test_clion ${SOURCE_FILES})

add_executable(Keldysh_mfRG
        main.cpp
        state.h
        vertex.h
        selfenergy.h
        propagator.h
        integrator.h
        data_structures.h
        frequency_grid.h
        parameters.h
        util.h
        OldFiles/save_hdf5.h
        Keldysh_symmetries.h
        OldFiles/susceptibility.h
        free_propagator_functions.h
        OldFiles/a_bubble.h
        OldFiles/p_bubble.h
        OldFiles/t_bubble.h loop.h a_vertex.h p_vertex.h t_vertex.h
        include/domain.hpp
        include/paid.hpp
        include/queue.hpp
        include/rule.hpp internal_symmetries.h
        bubbles.h diagrammatic_combinations.h)

include_directories(/usr/include/mpi/)

link_directories(/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib)

target_link_libraries(Keldysh_mfRG hdf5)
target_link_libraries(Keldysh_mfRG hdf5_cpp)
