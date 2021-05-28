cmake_minimum_required(VERSION 3.10)
project(Keldysh_mfRG)
project(Keldysh_mfRG_testing)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

set(CMAKE_CXX_COMPILER mpicxx)

# added -fopenmp
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS} -fopenmp") # unnecessary??
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage") // for code coverage analysis (slows down the code a lot!)

#set(SOURCE_FILES main.cpp) # diagrammatic_combinations.h interpolations.h testFunctions.h correctionFunctions.h)
#add_executable(openmp_test_clion ${SOURCE_FILES})

# standard executable for integration tests and production runs
add_executable(Keldysh_mfRG
        main.cpp
        bethe-salpeter.h
        bubbles.h
        correctionFunctions.h
        data_structures.h
        flow.h
        fourier_trafo.h
        frequency_grid.h
        hdf5_routines.h
        integrator.h
        Integrator_NR/integrator_NR.h
        internal_symmetries.h
        Keldysh_symmetries.h
        KramersKronig.h
        loop.h
        mpi_setup.h
        parameters.h
        parquet_checks.h
        perturbation_theory.h
        postprocessing.h
        propagator.h
        right_hand_sides.h
        r_vertex.h
        selfenergy.h
        solvers.h
        state.h
        symmetry_transformations.h
        symmetry_table.h
        testFunctions.h
        util.h
        vertex.h
        write_data2file.h momentum_grid.h)

# executable for unit tests using Catch
add_executable(Keldysh_mfRG_testing
        unit_tests/unit_tests.cpp
        unit_tests/catch.hpp
        unit_tests/test_data_structures.h
        unit_tests/test_symmetry_transformations.h
        unit_tests/test_integrator.h
        unit_tests/test_PrecalculateBubble.h
        bubbles.h
        bethe-salpeter.h
        correctionFunctions.h
        data_structures.h
        flow.h
        fourier_trafo.h
        frequency_grid.h
        hdf5_routines.h
        integrator.h
        Integrator_NR/integrator_NR.h
        internal_symmetries.h
        Keldysh_symmetries.h
        KramersKronig.h
        loop.h
        mpi_setup.h
        parameters.h
        parquet_checks.h
        perturbation_theory.h
        postprocessing.h
        propagator.h
        right_hand_sides.h
        r_vertex.h
        selfenergy.h
        solvers.h
        state.h
        symmetry_transformations.h
        symmetry_table.h
        testFunctions.h
        util.h
        vertex.h
        write_data2file.h momentum_grid.h unit_tests/test_momentum_grid.h)

include_directories(/usr/include/mpi/)
include_directories(/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/include)

target_link_libraries(Keldysh_mfRG -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
target_link_libraries(Keldysh_mfRG -lgsl -lgslcblas)
target_link_libraries(Keldysh_mfRG fftw3)

target_link_libraries(Keldysh_mfRG_testing -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
target_link_libraries(Keldysh_mfRG_testing -lgsl -lgslcblas)
target_link_libraries(Keldysh_mfRG_testing fftw3)