#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <set>
#include <thread>

#include "integrand.hpp"
#include "paid.hpp"

#define M_PI 3.14159265358979323846


class functor {
 public:
  double a;
  functor(double a_in):a(a_in){};

  std::complex<double> operator()(double x) {
   	std::cout<<"a="<<a<<std::endl;
    return std::complex<double>(std::exp(-(x * x))*a, std::exp(-(x * x)));
  }
};


int main() {
  double end = 2;
  Domain1D<std::complex<double>> D(0, end);
  double reference_value = std::erf(end) / 2 * sqrt(M_PI);
  functor f(1.0);

  std::vector<PAIDInput> inputs;
  inputs.push_back( PAIDInput( D, f, 1 ) );

  PAID p(inputs);

  auto result = p.solve();
  std::cout << "result: " << std::real(result[1]) << " ref: " << reference_value << "\n";
}
