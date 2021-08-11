cmake_minimum_required(VERSION 3.10)
project(Keldysh_mfRG)
project(Keldysh_mfRG_testing)
project(FPP_mfRG)

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
        bubbles.h
        postprocessing/causality_FDT_checks.h
        correctionFunctions.h
        data_structures.h
        flow.h
        grids/frequency_grid.h
        utilities/hdf5_routines.h
        integrator/integrator.h
        integrator/integrator_NR.h
        symmetries/internal_symmetries.h
        symmetries/Keldysh_symmetries.h
        postprocessing/KramersKronig.h
        loop.h
        utilities/mpi_setup.h
        parameters.h
        parquet_solver.h
        perturbation_theory.h
        postprocessing/postprocessing.h
        propagator.h
        right_hand_sides.h
        r_vertex.h
        selfenergy.h
        ODE_solvers.h
        state.h
        symmetries/symmetry_transformations.h
        symmetries/symmetry_table.h
        tests/test_perturbation_theory.h
        utilities/util.h
        vertex.h
        utilities/write_data2file.h grids/momentum_grid.h tests/test_ODE.h integrator/old_integrators.h)

# executable for unit tests using Catch
add_executable(Keldysh_mfRG_testing
        tests/unit_tests/unit_tests.cpp
        tests/unit_tests/catch.hpp
        tests/unit_tests/test_data_structures.h
        tests/unit_tests/test_symmetry_transformations.h
        tests/unit_tests/test_integrator.h
        tests/test_PrecalculateBubble.h
        bubbles.h
        postprocessing/causality_FDT_checks.h
        correctionFunctions.h
        data_structures.h
        flow.h
        grids/frequency_grid.h
        utilities/hdf5_routines.h
        integrator/integrator.h
        integrator/integrator_NR.h
        symmetries/internal_symmetries.h
        symmetries/Keldysh_symmetries.h
        postprocessing/KramersKronig.h
        loop.h
        utilities/mpi_setup.h
        parameters.h
        parquet_solver.h
        perturbation_theory.h
        postprocessing/postprocessing.h
        propagator.h
        right_hand_sides.h
        r_vertex.h
        selfenergy.h
        ODE_solvers.h
        state.h
        symmetries/symmetry_transformations.h
        symmetries/symmetry_table.h
        tests/test_perturbation_theory.h
        utilities/util.h
        vertex.h
        utilities/write_data2file.h
        grids/momentum_grid.h
        tests/unit_tests/test_momentum_grid.h
        tests/test_ODE.h
        integrator/old_integrators.h)

# standard executable for integration tests and production runs
add_executable(FPP_mfRG
        FPP/FPP_main.cpp
        FPP/Momentum-integral-Bubble.h
        FPP/Coordinates.h
        FPP/Differential_Equation.h
        FPP/1D-integrals.h
        FPP/Monte-Carlo_Trial.h
        FPP/Ladder-approximation.h
        bubbles.h
        postprocessing/causality_FDT_checks.h
        correctionFunctions.h
        data_structures.h
        flow.h
        grids/frequency_grid.h
        utilities/hdf5_routines.h
        integrator/integrator.h
        integrator/integrator_NR.h
        symmetries/internal_symmetries.h
        symmetries/Keldysh_symmetries.h
        postprocessing/KramersKronig.h
        loop.h
        utilities/mpi_setup.h
        parameters.h
        parquet_solver.h
        perturbation_theory.h
        postprocessing/postprocessing.h
        propagator.h
        right_hand_sides.h
        r_vertex.h
        selfenergy.h
        ODE_solvers.h
        state.h
        symmetries/symmetry_transformations.h
        symmetries/symmetry_table.h
        tests/test_perturbation_theory.h
        utilities/util.h
        vertex.h
        utilities/write_data2file.h
        grids/momentum_grid.h
        grids/flow_grid.h
        tests/test_ODE.h
        integrator/old_integrators.h)

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