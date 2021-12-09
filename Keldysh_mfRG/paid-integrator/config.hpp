// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:
# ifndef  PAID_INTEGRATOR_CONFIG_HPP
# define PAID_INTEGRATOR_CONFIG_HPP

#pragma once

#include "types.hpp"

namespace paid {

struct PAIDConfig {
 public:
    PAIDConfig()
            : max_f_evals(1e+9), // 1e+9
              max_iterations(-1), // -1
              max_error(1e-10), // 1e-10
              relative_error(0), // 0
              integration_rule(IntegrationRule::ClenshawCurtis),
              order(8), // 64
              ntasks_per_iteration(25),   // 25
              check_every_iteration(200), // 200
              check_below_iteration(10),  // 10
              min_size(1e-14), // 1e-14
              min_error(1e-20),  // 1e-20
              correct_error(true),    // true
              keep_small(true) {} // true

    PAIDConfig(std::ptrdiff_t max_f_evals_, double max_error_, bool relative_error_, std::size_t order_)
            : max_f_evals(max_f_evals_), // 1e+9
              max_iterations(-1), // -1
              max_error(max_error_), // 1e-10
              relative_error(relative_error_), // 0
              integration_rule(IntegrationRule::ClenshawCurtis),
              order(order_), // 64
              ntasks_per_iteration(25),   // 25
              check_every_iteration(200), // 200
              check_below_iteration(10),  // 10
              min_size(1e-14), // 1e-14
              min_error(1e-20),  // 1e-20
              correct_error(true),    // true
              keep_small(true) {} // true

  // maxima or goals
  std::ptrdiff_t max_f_evals;
  std::ptrdiff_t max_iterations;
  double max_error;
  //std::function<bool(double)> check_error = [](double error) -> bool { return error < max_error; };
  bool relative_error;

  // numerical method config
  IntegrationRule integration_rule;
  std::size_t order;

  // algorithmic details
  std::size_t ntasks_per_iteration;
  std::ptrdiff_t check_every_iteration;
  std::ptrdiff_t check_below_iteration;
  double min_size;
  double min_error;
  bool correct_error;
  bool keep_small;  // TODO
};

}  // namespace paid

#endif