/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// This file is a part of PAID.
// Copyright (c) 2016-2018, Simulation Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany
// and
// Copyright (c) 2016-2018, Aachen Institute for Advanced Study in Computational
//   Engineering Science, RWTH Aachen University, Germany All rights reserved.
// License is TODO

#include <cmath>
#include <complex>
#include <functional>
#include <vector>

// A construct used for templating that 'unwraps' the complex type to its
// value_type, and does not change other types.
// What we want is that
// Base< std::complex< double > > -> double
// Base<               double   > -> double
template <class Q>
struct Base_Class {
  typedef Q type;
};

template <class Q>
struct Base_Class<std::complex<Q>> {
  typedef Q type;
};

template <typename Q>
using Base = typename Base_Class<Q>::type;

template <typename T>
struct IntegrationResult {
  std::size_t fevals;
  T value;
  Base<T> error;
};

#include "domain.hpp"
#include "queue.hpp"
#include "rule.hpp"
//
// typedef std::function<T(Base<T>, Base<T>)> F;

class PAIDInput {
 public:
  PAIDInput(Domain1D<std::complex<double>> d_,
            std::function<std::complex<double>(double)> f_, std::size_t idx_)
      : d(d_), f(f_), idx(idx_) {}
  Domain1D<std::complex<double>> d;
  std::function<std::complex<double>(double)> f;
  std::size_t idx;
};

class PAID {
 public:
  using T = Domain1D<std::complex<double>>;

  PAID(std::vector<PAIDInput> inputs) : fevals(0), rule(16) {
    for (auto input : inputs) {
      // construct task:
      auto new_f = input.d.transform(input.f);
      auto res = rule.apply(new_f);
      fevals += res.fevals;
      q.put({input.d, input.f, res.value, res.error, input.idx});
    }
  }

  std::map< std::size_t, T::value_type > solve() {
    int done = 0;

    while (!done && q.sum().first > 1e-12) {
      // if (i == 0) std::cout << "current error: " << q.sum().first << "\n";
      if (q.empty()) {
        done++;
        break;
      }
      auto task = q.pop();
      auto doms = task.d.split();

      for (auto dom : doms) {
        auto new_f = dom.transform(task.f);
        auto res = rule.apply(new_f);

        fevals += res.fevals;
        q.put({dom, task.f, res.value, res.error, task.idx});
      }
    }

    auto res = q.sum();
    return res.second;
  }

 private:
  Queue<T> q;
  ClenshawCurtis<typename T::value_type> rule;
  std::size_t fevals;
};
