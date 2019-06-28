/*
 * Useful functions such as measuring runtime
 */

#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>  // system time
#include <unistd.h>    // time delay
#include <string>

using namespace std;

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
  cout << "time elapsed: ";
  printf("%.3f", t-t0);
  cout << "s" << endl;
}

/* display time difference in seconds w.r.t. reference time, with microsecond precision */
void get_time(double t0, string prec) {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  if (prec == "us") {
    long int us = tp.tv_sec * 1000000 + tp.tv_usec;
    double t = us / 1000000.;
    cout << "time elapsed: ";
    printf("%.6f", t-t0);
    cout << "s" << endl;
  }
  else {
    long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
    double t = ms / 1000.;
    cout << "time elapsed: ";
    printf("%.3f", t-t0);
    cout << "s" << endl;
  }
}

#endif // UTIL_H
