cmake_minimum_required(VERSION 3.14)
project(Keldysh_SelfEnergy_SOPT_2_0)

set(CMAKE_CXX_STANDARD 14)

add_executable(Keldysh_SelfEnergy_SOPT_2_0 main.cpp Propagator.h util.h)