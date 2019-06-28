cmake_minimum_required(VERSION 3.14)
project(Integration_for_Comparing__Debugging_)

set(CMAKE_CXX_STANDARD 14)

include_directories(.
        /include
        /tests)

add_executable(Integration_for_Comparing__Debugging_
        approxxpp.cpp
        approxxpp.h
        blnla.cpp
        blnla.h
        CMakeLists.txt
        integrate_new.h
        main.cpp
        matrix.cpp
        matrix.h
        patterson_set.cpp
        patterson_set.h
        util.h
        tests/integrand.hpp
        tests/test.cpp
        include/domain.hpp
        include/paid.hpp
        include/queue.hpp
        include/rule.hpp)
