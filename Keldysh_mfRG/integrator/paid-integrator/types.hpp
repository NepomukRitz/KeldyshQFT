// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:
# ifndef  PAID_INTEGRATOR_TYPES_HPP
# define PAID_INTEGRATOR_TYPES_HPP
#pragma once

#include <cstdio>

namespace paid {

// error_type depends on the compiler
#define PAID_USE_QUAD 1

#if defined(PAID_USE_QUAD) && defined(__INTEL_COMPILER)
using error_type = _Quad;
#elif defined(PAID_USE_QUAD) && (defined(__GNUC__) || defined(__clang__))
using error_type = __float128;
#else
using error_type = double;
#endif

enum IntegrationRule { ClenshawCurtis };

enum TerminationReason {
  Tolerance_Reached,
  No_Work_Left,
  Max_Evals_Reached,
  Max_Iter_Reached,
  Unknown,
  PAID_Not_Valid
};

template <typename T>
struct IntegrationResult {
  std::size_t fevals_;
  T value;
  error_type error;
};

// substitution of a general variable to x running in the interval of integration rule
struct AffineTransform1D {
public:
  double inline operator()(double x) const noexcept { return x * dmc + dpc; }

  double dmc;
  double dpc;
};

// transforms an n-dimensional index i = (i_1, ..., i_n) into one index
// i = N^(n-1)*i_1 + N^(n-2)*i_2 + ... + i_n
template<std::size_t Dim>
static std::size_t get_composite_index(std::size_t Nmax, std::array<std::size_t,Dim> i_vector) {
        std::size_t result = 0;
        std::size_t factor_j = 1;
        for (std::size_t j = 0; j < Dim; ++j){
            result += factor_j*i_vector[Dim-j-1];
            //std::cout << "j = " << j << ", i[" << Dim-j-1 << "] = " << i_vector[Dim-j-1] << ", factor_j = " << factor_j << ", i = " << result << "\n";
            factor_j *= Nmax;
        }
        return result;
};

// i_1 = i/N^(n-1), i_2 = (i-N^(n-1)*i_1)/N^(n-2), i_3 = (i-N^(n-1)*i_1-N^(n-2)*i_2)/N^(n-3), ..., i_n = (i - sum_{j=1}^{n-1}N^(n-j)i_j)/1
template<std::size_t Dim>
static std::array<std::size_t,Dim> get_each_index(std::size_t Nmax, std::size_t composite_index) {
        std::array<std::size_t,Dim> result;
        std::size_t factor_j = pow(Nmax,Dim-1);
        std::size_t i_j = 0;
        std::size_t i_subtracted = 0;

        for (std::size_t j = 0; j < Dim; ++j){
            i_j = (composite_index-i_subtracted)/factor_j;
            result[j] = i_j;
            i_subtracted += factor_j*i_j;
            factor_j /= Nmax;
        }
        return result;
};

}  // namespace paid
#endif