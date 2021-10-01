/*
 * Useful functions such as measuring runtime
 */

#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>  // system time
#include <chrono>      // system time
#include <unistd.h>    // time delay
#include <iostream>    // text input/output
#include "../data_structures.h"

#include "mpi_setup.h"


// print a time stamp in the following format: YYYY-MM-DD | hh-mm-ss |
void print_time_stamp() {
    time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()); // get time stamp
    tm loc = *localtime(&tt);                                            // convert time stamp to readable format

    // get elements of time stamp: Y, M, D, h, m, s
    int tms [6] {loc.tm_year + 1900, loc.tm_mon + 1, loc.tm_mday, loc.tm_hour, loc.tm_min, loc.tm_sec};
    // separators for readable time stamp format
    std::string separators [6] {"-", "-", " | ", ":", ":", " | "};

    for (int i=0; i<6; ++i) {
        if (tms[i] < 10) std::cout << "0";     // make sure to use two-digit format: add zero if necessary
        std::cout << tms[i] << separators[i];  // print time organized by separators
    }
}

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
double get_time() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms/1000.;
    return t;
}

// display time difference in seconds w.r.t. reference time, with millisecond precision
void get_time(double t0) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms / 1000.;
    if (MPI_FLAG) {
        if (mpi_world_rank() == 0) {
            std::cout << "time elapsed: ";
            printf("%.3f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
    else {
        std::cout << "time elapsed: ";
        printf("%.3f", t-t0);
        std::cout << "s" << std::endl;
    }
}

// display time difference in seconds w.r.t. reference time, with microsecond precision
void get_time(double t0, std::string prec) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    if (prec == "us") {
        long int us = tp.tv_sec * 1000000 + tp.tv_usec;
        double t = us / 1000000.;
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << "time elapsed: ";
                printf("%.6f", t - t0);
                std::cout << "s" << std::endl;
            }
        }
        else {
            std::cout << "time elapsed: ";
            printf("%.6f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
    else {
        long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        double t = ms / 1000.;
        if (MPI_FLAG) {
            if (mpi_world_rank() == 0) {
                std::cout << "time elapsed: ";
                printf("%.3f", t - t0);
                std::cout << "s" << std::endl;
            }
    }
        else {
            std::cout << "time elapsed: ";
            printf("%.3f", t - t0);
            std::cout << "s" << std::endl;
        }
    }
}



template<int degreeplus, typename Q>
inline auto lagrangePoly(const Q x, const double (&xs)[degreeplus], const Q (&ys) [degreeplus]) -> Q {
    Q result = 0.;

    double denominator, numerator;
    for (int i = 0; i < degreeplus; i++) {
        numerator = 1.;
        denominator = 1.;
        for (int k = 0; k < i; k++) {
            denominator *= (xs[i] - xs[k]);
            numerator *= (x - xs[k]);
        }
        for (int k = i+1; k < degreeplus; k++) {
            denominator *= (xs[i] - xs[k]);
            numerator *= (x - xs[k]);
        }

        result += ys[i] * numerator/denominator;
    }

    return result;
}


#endif // UTIL_H