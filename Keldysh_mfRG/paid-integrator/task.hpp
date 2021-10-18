// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <cstdio>
#include <iostream>

#include "domain.hpp"

namespace paid {
template <std::size_t N, typename F, typename T>
struct Task {
  using function_type = F;
  using value_type = T;
  static constexpr std::size_t dimension = N;

  Task() : d(0, 1){};  // we cheat a little bit

  explicit Task(Domain<N> d_, F f_, T val_, error_type error_,
                std::size_t idx_) noexcept
      : d(d_), f(f_), val(val_), error(error_ * d_.size()), idx(idx_) {
    assert(error >= 0);
  }

  Task(const Task<N, F, T>&) = default;               // Copy Constructor
  Task(Task<N, F, T>&&) = default;                    // Move Constructor
  Task& operator=(const Task<N, F, T>&) & = default;  // Copy assignement
  Task& operator=(Task<N, F, T>&&) & = default;       // Move assignment

  // is updated error smaller?
  bool operator<(const Task& rhs) const noexcept { return error < rhs.error; }

  Domain<N> d;
  F f;
  T val;
  error_type error;
  std::size_t idx;
};

// print “Task” as “[domain, value, error]”
template <std::size_t N, typename F, typename T>
std::ostream& operator<<(std::ostream& oss_, const Task<N, F, T>& rhs) {
  std::ostringstream oss;
  oss << "([" << rhs.d << "] \t" << rhs.val << "\t" << rhs.error << ")";
  oss_ << oss.str();
  return oss_;
}

}  // namespace paid
