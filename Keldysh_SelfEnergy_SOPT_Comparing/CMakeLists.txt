cmake_minimum_required(VERSION 3.14)
project(Keldysh_SelfEnergy_SOPT_Comparing)

set(CMAKE_CXX_STANDARD 14)

add_executable(Keldysh_SelfEnergy_SOPT_Comparing main.cpp Propagator.h util.h approxxpp.h approxxpp.cpp blnla.cpp blnla.h integrate_new.h matrix.cpp matrix.h patterson_set.cpp patterson_set.h)