cmake_minimum_required(VERSION 3.14)
project(Keldysh_mfRG)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

include_directories("$ gsl-2.5/")

add_executable(Keldysh_mfRG
        main.cpp
        state.h
        vertex.h
        selfenergy.h
        propagator.h
        bubbles.h
        integrator.h
        data_structures.h
        frequency_grid.h
        parameters.h
        util.h
        save_hdf5.h
        Keldysh_symmetries.h
        susceptibility.h
        free_propagator_functions.h
        a_bubble.h
        p_bubble.h
        t_bubble.h loop.h a_vertex.h p_vertex.h t_vertex.h
        include/domain.hpp
        include/paid.hpp
        include/queue.hpp
        include/rule.hpp internal_symmetries.h)