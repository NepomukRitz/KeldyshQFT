/*
 * Useful functions such as measuring runtime
 */

#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>  // system time
#include <unistd.h>    // time delay
#include <string>
#include <iostream>
#include "mpi_setup.h"

using namespace std;

/* print string in standard output */
void print(string s) {
    if (mpi_world_rank() == 0) {
        cout << s;
    }
}

/* print char in standard output */
void print(char c) {
    if (mpi_world_rank() == 0) {
        cout << c;
    }
}

/* print double in standard output */
void print(double d) {
    if (mpi_world_rank() == 0) {
        cout << d;
    }
}

/* print string in standard output and add new line */
void print(string s, bool endline) {
    if (mpi_world_rank() == 0) {
        cout << s;
        if (endline) cout << endl;
    }
}

/* return time stamp in seconds with millisecond precision */
double get_time() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms/1000.;
    return t;
}

/* display time difference in seconds w.r.t. reference time, with millisecond precision */
void get_time(double t0) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms/1000.;
    print("time elapsed: ");
    if (mpi_world_rank() == 0) printf("%.3f", t-t0);
    print("s", true);
}

/* display time difference in seconds w.r.t. reference time, with microsecond precision */
void get_time(double t0, std::string prec) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    if (prec == "us") {
        long int us = tp.tv_sec * 1000000 + tp.tv_usec;
        double t = us / 1000000.;
        print("time elapsed: ");
        if (mpi_world_rank() == 0) printf("%.6f", t-t0);
        print("s", true);
    }
    else {
        long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        double t = ms / 1000.;
        cout << "time elapsed: ";
        if (mpi_world_rank() == 0) printf("%.3f", t-t0);
        print("s", true);
    }
}

#endif // UTIL_H