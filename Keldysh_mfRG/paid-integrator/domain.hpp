// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <assert.h>
#include <array>
#include <ostream>
#include <sstream>

#include "types.hpp"

/*
template <typename T>
class Domain {
 private:
  virtual std::vector<Domain> split();
  virtual Base<T> size();
  virtual ~Domain();
};
*/

/*
template <std::size_t N>
class Domain {
 public:
  using T = std::array<std::pair<double, double>, N>;

  Domain(T arr) : dims_(arr) {}

  // Domain(T::iterator begin, T::iterator end) {
  //   std::copy(begin, end, dims_.begin());
  // }

  Domain(Domain<N - 1> dold, Domain<1> dnew) {
    std::copy(dold.cbegin(), dold.cend(), dims_.begin());
    dims_.back() = dnew.dims_[0];
  }

  std::vector<Domain<N> > split(Domain<N> dom) {
    std::array<std::pair<double, double>,N-1> oldarray;
    std::copy(dims_.cbegin(), dims_.cend(), oldarray.begin());

    Domain<N - 1> dold(oldarray);


    Domain<1> dnew{dims_.back()};
    auto dolds = dold.split();
    std::vector<Domain<N> > results;
    for (auto item : dolds) {
      results.push_back(Domain<N>(item, dnew));
    }
    return results;
  }

 private:
  T dims_;
};
*/

namespace paid {

template <std::size_t N>
class Domain {};

template <>
class Domain<1> {
 public:
  Domain(double left, double right) : left_(left), right_(right) {
    assert(left < right);
  }

  Domain(const Domain&) = default;               // Copy Constructor
  Domain(Domain&&) = default;                    // Move Constructor
  Domain& operator=(const Domain&) & = default;  // Copy Assigment
  Domain& operator=(Domain&&) & = default;       // Move Assigment

  std::array<Domain<1>, 2> split() const {
    auto length = right_ - left_;
    auto middle = left_ + length / 2;
    return std::array<Domain<1>, 2>{
        {Domain{left_, middle}, Domain{middle, right_}}};
  }

  // transforms general integration interval to specific one from the integration rule
  std::array<AffineTransform1D, 1> getTransform() const {
    double dmc = 0.5 * (right_ - left_);
    double dpc = 0.5 * (right_ + left_);
    return std::array<AffineTransform1D, 1>{{AffineTransform1D{dmc, dpc}}};
  }

  // width of the integration interval
  double size() const {
    assert(right_ - left_ > 0);
    return right_ - left_;
  }

  double left_;
  double right_;
};

// print "Domain" as "(left,right)"
std::ostream& operator<<(std::ostream& oss_, const Domain<1>& rhs) {
  std::ostringstream oss;
  oss << rhs.left_ << "," << rhs.right_ << " (" << rhs.size() << ")";
  oss_ << oss.str();
  return oss_;
}

}  // namespace paid
