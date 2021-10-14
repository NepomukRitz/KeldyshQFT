// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <assert.h>
#include <cmath>
#include <vector>

#include "invoker.hpp"
#include "types.hpp"

namespace paid {

// an integration rule
/*
class IntegrationRule {
 public:
};
*/

class ClenshawCurtis {
 public:
  explicit ClenshawCurtis(std::size_t N) noexcept
      : N_(N), w(N + 1), w2(2 * N + 1), t(N + 1), t2(2 * N + 1) {
    std::size_t j, k;

    std::vector<std::vector<double>> z(N, std::vector<double>(N + 1));
    std::vector<std::vector<double>> z2(2 * N, std::vector<double>(2 * N + 1));

    // z, z2, just notation to help to compute the weights
    for (j = 0; j < N_; ++j)
      for (k = 0; k <= N_; ++k) {
        z[j][k] = std::cos(k * j * M_PI / N_) / N_;
        if (k != 0 && k != N_) z[j][k] *= 2.0;
        if (j == 0)
          z[j][k] /= 2.0;
        else if (j % 2 == 0)
          z[j][k] /= ((j + 1) * (1 - j));
      }

    for (j = 0; j < 2 * N_; ++j)
      for (k = 0; k <= 2 * N_; ++k) {
        z2[j][k] = std::cos(k * j * M_PI / (2 * N_)) / (2 * N_);
        if (k != 0 && k != (2 * N_)) z2[j][k] *= 2.0;
        if (j == 0)
          z2[j][k] /= 2.0;
        else if (j % 2 == 0)
          z2[j][k] /= ((j + 1) * (1 - j));
      }

    // weights for N_ and 2N_
    for (k = 0; k <= N_; ++k) {
      w[k] = 0.0;
      for (j = 0; j < N_; j += 2) w[k] += z[j][k];
    }

    for (k = 0; k <= 2 * N_; ++k) {
      w2[k] = 0.0;
      for (j = 0; j < 2 * N_; j += 2) w2[k] += z2[j][k];
    }

    std::size_t npts = N_;

    for (k = 0; k <= npts; ++k) {
      t[k] = std::cos((M_PI * k) / static_cast<double>(npts));
      if (k == 0 || k == npts) {
        w[k] = 1.0 / (npts * npts - 1.0);
      } else {
        w[k] = 1.0 + std::cos(k * M_PI) / (1.0 - npts * npts);
        for (j = 1; j <= npts / 2 - 1; j++) {
          w[k] += (2.0 / (1.0 - 4.0 * j * j)) *
                  std::cos((2.0 * j) * (M_PI * k) / npts);
        }
        w[k] *= 2.0 / npts;
      }
    }

    npts = 2 * N_;

    for (k = 0; k <= npts; ++k) {
      t2[k] = std::cos((M_PI * k) / static_cast<double>(npts));
      if (k == 0 || k == npts) {
        w2[k] = 1.0 / (npts * npts - 1.0);
      } else {
        w2[k] = 1.0 + std::cos(k * M_PI) / (1.0 - npts * npts);
        for (j = 1; j <= npts / 2 - 1; j++) {
          w2[k] += (2.0 / (1.0 - 4.0 * j * j)) *
                   std::cos((2.0 * j) * (M_PI * k) / npts);
        }
        w2[k] *= 2.0 / npts;
      }
    }
  }

  template <std::size_t N, typename F, typename T, typename... Args>
  IntegrationResult<T> apply(const F& f,
                             const std::array<AffineTransform1D, N>& af) const
      noexcept {
    T val1;
    T val2;
    const std::size_t nEvals = 2 * (N_ + 1);
    std::vector<T> cache_(N_ + 1);  // TODO

    for (std::size_t k = 0; k <= N_; ++k) {
      auto node = af[0](t[k]);
      cache_[k] =
          Invoker<N, F, T, Args...>::apply(f, std::array<double, 1>{{node}});
      val1 += w[k] * cache_[k];
    }

    val2 = 0;
    for (std::size_t k = 0; k <= 2 * N_; ++k) {
      if (k % 2 == 0) {
        val2 += w2[k] * cache_[k / 2];
        assert(t[k / 2] == t2[k]);
      } else {
        // val2 += w2[k] * f(af(t2[k]));
        auto node = af[0](t2[k]);
        val2 += w2[k] * Invoker<N, F, T, Args...>::apply(
                            f, std::array<double, 1>{{node}});
      }
    }

    assert(!std::isnan(std::abs(val2)));
    assert(!std::isnan(std::abs(val1)));

    return {nEvals, val2 * (af[0](1) - af[0](0)), std::abs(val2 - val1)};
  };

 private:
  std::size_t N_;
  std::vector<double> w;
  std::vector<double> w2;
  std::vector<double> t;
  std::vector<double> t2;
};

}  // namespace paid
