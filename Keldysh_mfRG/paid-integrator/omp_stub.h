// This file is a part of PAID.
// Copyright (c) 2015-2021, Simulation and Data Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany. All rights reserved.
// License is 3-clause BSD:

#pragma once

// A stub for omp
// Markus insists that this a good idea.

int omp_get_num_threads() { return 1; }
int omp_get_max_threads() { return 1; }
int omp_get_thread_num() { return 0; }
