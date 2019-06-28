cmake_minimum_required(VERSION 3.14)
project(one_loop_fRG_wPAID)

set(CMAKE_CXX_STANDARD 14)

include_directories(.
        /include)

add_executable(one_loop_fRG_wPAID main.cpp
        include/domain.hpp
        include/paid.hpp
        include/queue.hpp
        include/rule.hpp
        util.h)