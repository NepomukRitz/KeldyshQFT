// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:
# ifndef  PAID_INTEGRATOR_TASK_HPP
# define PAID_INTEGRATOR_TASK_HPP
#pragma once

#include <cstdio>
#include <iostream>

#include "domain.hpp"

namespace paid {
template <std::size_t Dim, typename F, typename T, typename Key>
struct Task {
  //using function_type = F;
  //using value_type = T;
  static constexpr std::size_t dimension = Dim;

  Task() : d({0}, {1}){};  // we cheat a little bit

  explicit Task(Domain<Dim> d_, F f_, T val_, error_type error_,
                Key idx_) noexcept
      : d(d_), f(f_), val(val_), error(error_ * d_.size()), idx(idx_) {
    assert(error >= 0);
  }

  Task(const Task<Dim, F, T, Key>&) = default;               // Copy Constructor
  Task(Task<Dim, F, T, Key>&&) = default;                    // Move Constructor
  //Task& operator=(const Task<Dim, F, T, Key>&) & = default;  // Copy assignement
  Task& operator=(Task<Dim, F, T, Key>&&) & = default;       // Move assignment
  Task& operator=(const Task<Dim, F, T, Key>& other){ // Copy assignment
        this->d = other.d;
        this->val = other.val;
        this->error = other.error;
        this->idx = other.idx;
        this->f = other.f;
        //do not copy f

        return *this;
  }

  // is updated error smaller?
  bool operator<(const Task& rhs) const noexcept { return error < rhs.error; }

  Domain<Dim> d;
  F f;
  T val;
  error_type error;
  Key idx;
};

// print “Task” as “[domain, value, error]”
template <std::size_t Dim, typename F, typename T, typename Key>
std::ostream& operator<<(std::ostream& oss_, const Task<Dim, F, T, Key>& rhs) {
  std::ostringstream oss;
  oss << "([" << rhs.d << "] \t" << rhs.val << "\t" << rhs.error << ")";
  oss_ << oss.str();
  return oss_;
}

}  // namespace paid
#endif