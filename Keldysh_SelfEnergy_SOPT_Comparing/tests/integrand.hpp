#pragma once

#include <cmath>
#include <complex>

// This files contains example integrands.

#define W 1e-6
#define LAMBDA 1e-6
#define Q1 0
#define Q2 0
#define BIGW 0

double d1 = 1.3;
double d2 = 2.75;
double m1 = 0.4;
double m2 = 1.8;

std::complex<double> g(double omega) {
  std::complex<double> i(0.0, 1.0);
  return i * sqrt(omega * omega + LAMBDA * LAMBDA);
}

std::complex<double> d_lambda_g(double omega) {
  std::complex<double> i(0.0, 1.0);
  return i * 2.0 * LAMBDA / sqrt(omega * omega + LAMBDA * LAMBDA);
}

double e(double k1, double k2) { return 2.0 * (cos(k1) + cos(k2)); }

double f1(double x) { return 1 - 3 * x / 25; }

double f2(double x) { return 1 + 3 * x / 25; }

double h1(double x) { return (-5) * ((4.5 - x) / (5.5 - x)); }

double h2(double x) { return 5 * ((4.5 + x) / (5.5 + x)); }

std::complex<double> Sigma(double omega, double eps) {
  std::complex<double> i(0.0, 1.0);
  return m1 / (eps + i * (omega + d1 / 2.0)) +
         m2 * f1(eps) / (h1(eps) + i * (omega + d2 / 2.0)) +
         m2 * f2(eps) / (h2(eps) + i * (omega + d2 / 2.0));
}

std::complex<double> G(double q1, double q2, double omega) {
  return 1.0 / (g(omega) - e(q1, q2) - Sigma(omega, e(q1, q2)));
}

std::complex<double> S(double q1, double q2, double omega) {
  return (-1.0) * G(q1, q2, omega) * G(q1, q2, omega) * d_lambda_g(omega);
}

// theta is 1 on whole domain, for testing
double theta(double q1, double q2) { return 1.0; }

/*
std::complex<double> integrand( double q1, double q2 )
{
    return S(q1,q2,W)*G(Q1-q1,Q2-q2,BIGW-W)*theta(q1,q2);
}
*/

//------------------------------------------------------------
std::complex<double> frac(double x, double a, double eps) {
  return std::complex<double>(1, 0) / (std::complex<double>(x + a, eps));
}

double realFrac(double x, double a, double eps) { return 1 / (x + a + eps); }

std::complex<double> evilFrac(double x, double y, double a, double eps) {
  return std::complex<double>(1, 0) / (std::complex<double>(x + y + a, eps));
}

double integrand(double x, double y) {
  // return 1/( std::exp( x*x + y*y ) - 0.2 );// * realFrac( y, -M_PI/3, 1e-1);
  return std::abs(std::sin(std::sin(x + y)));
}

std::complex<double> complex_integrand(double x, double y) {
  return frac(x, -M_PI / 2, 1e-2) * frac(y, -M_PI / 3, 1e-3);
}
