// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

#include <complex>
#include <functional>

#include "paid.hpp"

namespace lmu {

// using Key = std::size_t;
// using T_type = std::complex<double>;

using F_type = std::function<std::complex<double>(double)>;
static constexpr std::size_t Dim = 1;

template <typename T_type, typename Key>
using PAID = paid::PAID<Dim, F_type, T_type, Key, double>;

template <typename Key>
using PAIDInput = paid::PAIDInput<Dim, F_type, Key>;

using PAIDConfig = paid::PAIDConfig;

}  // namespace paid
