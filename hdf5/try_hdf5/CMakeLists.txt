cmake_minimum_required(VERSION 3.14)
project(try_hdf5)

set(CMAKE_CXX_STANDARD 11)

include_directories(/usr/local/Cellar/hdf5/1.10.5_1/include)

add_executable(try_hdf5 main.cpp save_hdf5.h)