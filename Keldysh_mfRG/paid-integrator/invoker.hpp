// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:
# ifndef  PAID_INTEGRATOR_INVOKER_HPP
# define PAID_INTEGRATOR_INVOKER_HPP

#pragma once

#include <cstdio>
#include <array>

namespace paid {

template <std::size_t Dim, typename F, typename T, typename... Args>
struct Invoker { // invoker of the function without an applicator
  static T apply(const F&, std::array<double, Dim>) {
    /// no specialization exists!.
    /// Probably your function dimension N is too high
    return F::this_type_is_missing_an_applicator();
  }
};

template <typename F, typename T>
struct Invoker<1, F, T, double> { // invoker of the integrand function
  static T apply(const F& f, std::array<double, 1> input) {
    return f(input[0]);
  }
};

template <std::size_t Dim, typename F, typename T>
struct Invoker<Dim, F, T, std::array<double, Dim>> { // invoker of higher-dimensional integrand function
  static T apply(const F& f, std::array<double, Dim> input) { return f(input); }
};

}  // namespace paid

#endif