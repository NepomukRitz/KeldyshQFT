// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

# ifndef  PAID_INTEGRATOR_DOMAIN_HPP
# define PAID_INTEGRATOR_DOMAIN_HPP

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
template <std::size_t Dim>
class Domain {
 public:
  using T = std::array<std::pair<double, double>, Dim>;

  Domain(T arr) : dims_(arr) {}

  // Domain(T::iterator begin, T::iterator end) {
  //   std::copy(begin, end, dims_.begin());
  // }

  Domain(Domain<Dim - 1> dold, Domain<1> dnew) {
    std::copy(dold.cbegin(), dold.cend(), dims_.begin());
    dims_.back() = dnew.dims_[0];
  }

  std::vector<Domain<Dim> > split(Domain<Dim> dom) {
    std::array<std::pair<double, double>,Dim-1> oldarray;
    std::copy(dims_.cbegin(), dims_.cend(), oldarray.begin());

    Domain<Dim - 1> dold(oldarray);


    Domain<1> dnew{dims_.back()};
    auto dolds = dold.split();
    std::vector<Domain<Dim> > results;
    for (auto item : dolds) {
      results.push_back(Domain<Dim>(item, dnew));
    }
    return results;
  }

 private:
  T dims_;
};
*/

namespace paid {
/*
template <std::size_t Dim>
class Domain {};

// 1D-case
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
*/

namespace {
    template <std::size_t n> constexpr std::size_t power_int(std::size_t x) {
        std::size_t result = 1;
        for (std::size_t n_counter = 0; n_counter < n; ++n_counter) {
            result *= x;
        }
        return result;
    }

    template <> constexpr std::size_t power_int<1>(std::size_t x) {
        return x;
    }
}

// general case
template <std::size_t Dim>
class Domain {
public:
    Domain(std::array<double,Dim> left, std::array<double,Dim> right) : left_(left), right_(right) {
        for(std::size_t n = 0; n < Dim; ++n){
            assert(left[n] < right[n]);
        }
    }

    Domain() {};                                   // Standard Constructor
    Domain(const Domain&) = default;               // Copy Constructor
    Domain(Domain&&) = default;                    // Move Constructor
    Domain& operator=(const Domain&) & = default;  // Copy Assigment
    Domain& operator=(Domain&&) & = default;       // Move Assigment

    // split interval of domain
    std::array<Domain<Dim>, power_int<Dim>(2)> split() const {
        std::array<double, Dim> length, middle;
        for (std::size_t n = 0; n < Dim; ++n) {
            length[n] = right_[n] - left_[n];
            middle[n] = left_[n] + length[n] / 2;
        }
        std::array<Domain<Dim>, power_int<Dim>(2)> split_domain;
        std::array<std::size_t, Dim> k_vector;
        std::array<double, Dim> left_helper, right_helper;
        for (std::size_t k = 0; k < power_int<Dim>(2); ++k) {
            left_helper = left_;
            right_helper = right_;
            k_vector = get_each_index<Dim>(2, k);
            for (std::size_t n = 0; n < Dim; ++n) {
                if (k_vector[n] == 0) {
                    left_helper[n] = left_[n];
                    right_helper[n] = middle[n];
                } else {
                    left_helper[n] = middle[n];
                    right_helper[n] = right_[n];
                }
            };
            split_domain[k] = Domain{left_helper, right_helper};
        };
        return split_domain;
    }

    // transforms general integration interval to specific one from the integration rule
    std::array<AffineTransform1D, Dim> getTransform() const {
        std::array<double, Dim> dmc, dpc;
        std::array<AffineTransform1D, Dim> AffineTransforms;
        for (std::size_t n = 0; n < Dim; ++n) {
            dmc[n] = 0.5 * (right_[n] - left_[n]);
            dpc[n] = 0.5 * (right_[n] + left_[n]);
            AffineTransforms[n] = AffineTransform1D{dmc[n], dpc[n]};
        }
        return AffineTransforms;
    }

    // width of the integration interval
    double size() const {
        double size_return = 1;
        for (std::size_t n = 0; n < Dim; ++n) {
            assert(right_[n] - left_[n] > 0);
            size_return *= (right_[n] - left_[n]);
        }
        return size_return;
    }

    std::array<double,Dim> left_;
    std::array<double,Dim> right_;
};

// 2D-case
/*
template <>
class Domain<2> {
public:
    Domain(double leftx, double rightx, double lefty, double righty) : leftx_(leftx), rightx_(rightx), lefty_(lefty), righty_(righty) {
        assert(leftx < rightx);
        assert(lefty < righty);
    }

    Domain(const Domain&) = default;               // Copy Constructor
    Domain(Domain&&) = default;                    // Move Constructor
    Domain& operator=(const Domain&) & = default;  // Copy Assigment
    Domain& operator=(Domain&&) & = default;       // Move Assigment

    std::array<Domain<2>, 4> split() const {
        auto lengthx = rightx_ - leftx_;
        auto middlex = leftx_ + lengthx / 2;
        auto lengthy = righty_ - lefty_;
        auto middley = lefty_ + lengthy / 2;
        return std::array<Domain<2>, 4>{
            {Domain{leftx_, middlex,lefty_,middley}, Domain{middlex, rightx_,lefty_,middley},
                    Domain{leftx_, middlex,middley,righty_}, Domain{middlex, rightx_,middley,righty_}}};
    }

    // transforms general integration interval to specific one from the integration rule
    std::array<AffineTransform1D, 2> getTransform() const {
        double dmcx = 0.5 * (rightx_ - leftx_);
        double dpcx = 0.5 * (rightx_ + leftx_);
        double dmcy = 0.5 * (righty_ - lefty_);
        double dpcy = 0.5 * (righty_ + lefty_);
        return std::array<AffineTransform1D, 2>{{AffineTransform1D{dmcx, dpcx},AffineTransform1D{dmcy, dpcy}}};
    }

    // width of the integration interval
    double size() const {
        assert(rightx_ - leftx_ > 0);
        assert(righty_ - lefty_ > 0);
        return (rightx_ - leftx_)*(righty_ - lefty_);
    }

    double leftx_;
    double rightx_;
    double lefty_;
    double righty_;
};
 */

// print "Domain<1>" as "left,right(size)"
template <std::size_t Dim>
std::ostream& operator<<(std::ostream& oss_, const Domain<Dim>& rhs) {
  std::ostringstream oss;
  for (std::size_t n = 0; n < Dim; ++n) {
      oss << "(" << rhs.left_[n] << "," << rhs.right_[n] << ")";
  }
  oss << " (" << rhs.size() << ")";
  oss_ << oss.str();
  return oss_;
}

/*
// print "Domain<2>" as "(leftx,rightx),(lefty,righty)(size)"
    std::ostream& operator<<(std::ostream& oss_, const Domain<2>& rhs) {
        std::ostringstream oss;
        oss << "(" << rhs.leftx_ << "," << rhs.rightx_ << "), " <<  "(" << rhs.lefty_ << "," << rhs.righty_ << "),(" << rhs.size() << ")";
        oss_ << oss.str();
        return oss_;
    }
*/

}  // namespace paid

#endif