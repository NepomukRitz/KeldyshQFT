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

#set(SOURCE_FILES main.cpp) # diagrammatic_combinations.h interpolations.h testFunctions.h correctionFunctions.h)
#add_executable(openmp_test_clion ${SOURCE_FILES})

add_executable(Keldysh_mfRG
        main.cpp
        a_vertex.h
        bubbles.h
        correctionFunctions.h
        data_structures.h
        diagrammatic_combinations.h
        flow.h
        fourier_trafo.h
        frequency_grid.h
        hdf5_routines.h
        integrator.h
        internal_symmetries.h
        Keldysh_symmetries.h
        loop.h
        mpi_setup.h
        parameters.h
        perturbation_theory.h
        propagator.h
        p_vertex.h
        right_hand_sides.h
        selfenergy.h
        solvers.h
        state.h
        testFunctions.h
        t_vertex.h
        util.h
        vertex.h
        write_data2file.h Integrator_NR/integrator_NR.h)

add_executable(Keldysh_mfRG_testing
        unit_tests/unit_tests.cpp
        a_vertex.h
        bubbles.h
        correctionFunctions.h
        data_structures.h
        unit_tests/test_data_structures.h
        diagrammatic_combinations.h
        flow.h
        fourier_trafo.h
        frequency_grid.h
        hdf5_routines.h
        integrator.h
        internal_symmetries.h
        Keldysh_symmetries.h
        loop.h
        mpi_setup.h
        parameters.h
        perturbation_theory.h
        propagator.h
        p_vertex.h
        right_hand_sides.h
        selfenergy.h
        solvers.h
        state.h
        testFunctions.h
        t_vertex.h
        util.h
        vertex.h
        write_data2file.h Integrator_NR/integrator_NR.h)

include_directories(/usr/include/mpi/)

target_link_libraries(Keldysh_mfRG -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
target_link_libraries(Keldysh_mfRG -lgsl -lgslcblas)
target_link_libraries(Keldysh_mfRG fftw3)

target_link_libraries(Keldysh_mfRG_testing -L/software/opt/bionic/x86_64/hdf5/1.10.5-gcc/lib -lhdf5 -lhdf5_cpp)
target_link_libraries(Keldysh_mfRG_testing -lgsl -lgslcblas)
target_link_libraries(Keldysh_mfRG_testing fftw3)