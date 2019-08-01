cmake_minimum_required(VERSION 3.14)
project(Keldysh_mfRG)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

add_executable(Keldysh_mfRG main.cpp vertex.h selfenergy.h data_structures.h parameters.h util.h)