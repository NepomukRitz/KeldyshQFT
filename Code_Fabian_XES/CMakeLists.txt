cmake_minimum_required(VERSION 3.10)
project(Code_Fabian_XES)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Code_Fabian_XES
        plots/.gitkeep
        plots/xes_plots.ipynb
        plots/xes_plots.pdf
        .gitkeep
        myfuncs.h
        myheader.h
        mymain.cpp
        mytests.h)
