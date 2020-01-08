//
// Created by E.Walter on 12/19/19.
//

#ifndef KELDYSH_MFRG_MPI_SETUP_H
#define KELDYSH_MFRG_MPI_SETUP_H

#include "data_structures.h"
#include "parameters.h"
//#include <omp.h>
#include <mpi.h>

// Get the rank(ID) of the current process
int mpi_world_rank() {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}

// Get the number of processes
int mpi_world_size() {
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}

template <typename Q>
vec<Q> mpi_initialize_buffer(int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    vec<Q> buffer (n_omp*(n_mpi/world_size+1));
    return buffer;
}

template <typename Q>
vec<Q> mpi_initialize_result(int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    vec<Q> result (n_omp * (n_mpi - (n_mpi % world_size) + world_size));
    return result;
}

void mpi_collect(vec<comp>& buffer, vec<comp>& result, int n_mpi, int n_omp) {
    int world_size = mpi_world_size();

    MPI_Allgather(&buffer[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_COMPLEX,
                  &result[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_COMPLEX, MPI_COMM_WORLD);
}

void mpi_collect(vec<double>& buffer, vec<double>& result, int n_mpi, int n_omp) {
    int world_size = mpi_world_size();

    MPI_Allgather(&buffer[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_DOUBLE,
                  &result[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
}

template <typename Q>
vec<Q> mpi_reorder_result(vec<Q>& result, int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    vec<Q> ordered_result;
    for (int i=0; i<n_mpi/world_size; ++i) {
        for (int j=0; j<world_size; ++j) {
            typename vec<Q>::const_iterator first = result.begin() + (j*(n_mpi/world_size+1) + i) * n_omp;
            typename vec<Q>::const_iterator last = first + n_omp;
            ordered_result.insert(ordered_result.end(), first, last);
        }
    }
    for (int j=0; j<world_size; ++j) {
        if (n_mpi % world_size > j) {
            typename vec<Q>::const_iterator first = result.begin() + (j*(n_mpi/world_size+1) + n_mpi/world_size) * n_omp;
            typename vec<Q>::const_iterator last = first + n_omp;
            ordered_result.insert(ordered_result.end(), first, last);
        }
    }
    return ordered_result;
}

#endif //KELDYSH_MFRG_MPI_SETUP_H
