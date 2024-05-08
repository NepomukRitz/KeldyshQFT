/*
 * Useful functions such as measuring runtime
 */

#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>  // system time
#include <chrono>      // system time
#include <unistd.h>    // time delay
#include <iostream>    // text input/output
#include <bits/stdc++.h>
#include <sys/stat.h>
#include "../data_structures.hpp"

#include "mpi_setup.hpp"

namespace utils {

    // print a time stamp in the following format: YYYY-MM-DD | hh-mm-ss |
    void print_time_stamp();

    // print any printable data in standard output and add new line
    template <typename T>
    void print(T s, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                print_time_stamp();
                std::cout << s;
                if (endline) std::cout << std::endl;
            }
        }
        else {
            print_time_stamp();
            std::cout << s;
            if (endline) std::cout << std::endl;
        }
    }

    // print any printable data in standard output and add new line (without time stamp)
    template <typename T>
    void print_add(T s, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << s;
                if (endline) std::cout << std::endl;
            }
        }
        else {
        std::cout << s;
        if (endline) std::cout << std::endl;
        }
    }

    // print two different data types in standard output and add new line
    template <typename T, typename U>
    void print(T t, U u, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                print_time_stamp();
                std::cout << t << u;
                if (endline) std::cout << std::endl;
            }
        }
        else {
            print_time_stamp();
            std::cout << t << u;
            if (endline) std::cout << std::endl;
        }
    }

    template<typename ...Args>
    void print(Args ... args) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                print_time_stamp();
                (std::cout << ... << args);
            }
        } else {
            print_time_stamp();
            (std::cout << ... << args);
        }
    }

    // print two different data types in standard output and add new line (without time stamp)
    template <typename T, typename U>
    void print_add(T t, U u, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << t << u;
                if (endline) std::cout << std::endl;
            }
        }
        else {
            std::cout << t << u;
            if (endline) std::cout << std::endl;
        }
    }

    // print three different data types in standard output and add new line
    template <typename T, typename U, typename V>
    void print(T t, U u, V v, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                print_time_stamp();
                std::cout << t << u << v;
                if (endline) std::cout << std::endl;
            }
        }
        else {
            print_time_stamp();
            std::cout << t << u << v;
            if (endline) std::cout << std::endl;
        }
    }

    // print three different data types in standard output and add new line (without time stamp)
    template <typename T, typename U, typename V>
    void print_add(T t, U u, V v, bool endline) {
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << t << u << v;
                if (endline) std::cout << std::endl;
            }
        }
        else {
            std::cout << t << u << v;
            if (endline) std::cout << std::endl;
        }
    }

    // print any printable data in standard output
    template <typename T>
    void print(T s) {
        print(s, false);
    }

    // print two different data types in standard output
    template <typename T, typename U>
    void print(T t, U u) {
        print(t, u, false);
    }

    // print three different data types in standard output
    template <typename T, typename U, typename V>
    void print(T t, U u, V v) {
        print(t, u, v, false);
    }

    // return time stamp in seconds with millisecond precision
    double get_time();

    // display time difference in seconds w.r.t. reference time, with millisecond precision
    void get_time(double t0);
    void get_cpu_hours(double t0);

    // display time difference in seconds w.r.t. reference time, with microsecond precision
    void get_time(double t0, std::string prec);


    void makedir(const std::string& dir_str);


    void check_input(const fRG_config& config);

    std::string generate_data_directory(std::string& job);

    std::string generate_filename(const fRG_config& config);

    void print_job_info(const fRG_config& config);

    void hello_world();

    //template <typename T, typename... Args>
    //constexpr std::array<T, sizeof...(Args)> to_array (const Args && ... args)
    //{ return {{ static_cast<T>(std::forward<Args>(args))... }}; }


    template <typename T, typename... Args>
    constexpr std::array<T, sizeof...(Args)> to_array (const Args & ... args)
    { return {{ static_cast<T>(args)... } }; }

    template<typename Head, typename ...Args>
    bool is_all_finite(Head head, Args ... args) {
        if constexpr (sizeof... (Args) == 0) {
            return std::isfinite(head);
        }
        else {
            return std::isfinite(head) and is_all_finite(args...);
        }
    }
}

#endif // UTIL_H