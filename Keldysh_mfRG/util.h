/*
 * Useful functions such as measuring runtime
 */

#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>  // system time
#include <unistd.h>    // time delay
#include <iostream>    // text input/output

#ifdef MPI_FLAG
#include "mpi_setup.h"
#endif

using namespace std;

// print any printable data in standard output and add new line
template <typename T>
void print(T s, bool endline) {
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0) {
        cout << s;
        if (endline) cout << endl;
    }
#else
    cout << s;
    if (endline) cout << endl;
#endif
}

// print two different data types in standard output and add new line
template <typename T, typename U>
void print(T t, U u, bool endline) {
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0) {
        cout << t << u;
        if (endline) cout << endl;
    }
#else
    cout << t << u;
    if (endline) cout << endl;
#endif
}

// print three different data types in standard output and add new line
template <typename T, typename U, typename V>
void print(T t, U u, V v, bool endline) {
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0) {
        cout << t << u << v;
        if (endline) cout << endl;
    }
#else
    cout << t << u << v;
    if (endline) cout << endl;
#endif
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
    double t = ms/1000.;
    print("time elapsed: ");
#ifdef MPI_FLAG
    if (mpi_world_rank() == 0) printf("%.3f", t-t0);
#else
    printf("%.3f", t-t0);
#endif
    print("s", true);
}

// display time difference in seconds w.r.t. reference time, with microsecond precision
void get_time(double t0, std::string prec) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    if (prec == "us") {
        long int us = tp.tv_sec * 1000000 + tp.tv_usec;
        double t = us / 1000000.;
        print("time elapsed: ");
#ifdef MPI_FLAG
        if (mpi_world_rank() == 0) printf("%.6f", t-t0);
#else
        printf("%.6f", t-t0);
#endif
        print("s", true);
    }
    else {
        long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        double t = ms / 1000.;
        cout << "time elapsed: ";
#ifdef MPI_FLAG
        if (mpi_world_rank() == 0) printf("%.3f", t-t0);
#else
        printf("%.3f", t-t0);
#endif
        print("s", true);
    }
}

#endif // UTIL_H