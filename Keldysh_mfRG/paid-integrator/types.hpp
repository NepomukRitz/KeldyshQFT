// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <cstdio>

namespace paid {

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
  double error;
};

struct AffineTransform1D {
 public:
  double inline operator()(double x) const noexcept { return x * dmc + dpc; }

  double dmc;
  double dpc;
};

}  // namespace paid
