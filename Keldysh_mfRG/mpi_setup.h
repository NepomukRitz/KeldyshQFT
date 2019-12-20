//
// Created by E.Walter on 12/19/19.
//

#ifndef KELDYSH_MFRG_MPI_SETUP_H
#define KELDYSH_MFRG_MPI_SETUP_H

#include "data_structures.h"
#include "parameters.h"
//#include <omp.h>
#include <mpi.h>

template <typename Q>
vec<Q> mpi_initialize_buffer(int n_para_mpi, int n_para_omp) {

    // Get the rank(ID) of the current process
    int world_rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    vec<Q> K3_buffer (n_para_omp*(n_para_mpi/world_size+1));

    return K3_buffer;

}

void mpi_collect(vec<comp>& K3_buffer, vec<comp>& K3_result, int n_para_mpi, int n_para_omp) {

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Allgather(&K3_buffer[0], static_cast<int>(n_para_omp*(n_para_mpi/world_size+1)), MPI_COMPLEX,
                  &K3_result[0], static_cast<int>(n_para_omp*(n_para_mpi/world_size+1)), MPI_COMPLEX, MPI_COMM_WORLD);
}

void mpi_collect(vec<double>& K3_buffer, vec<double>& K3_result, int n_para_mpi, int n_para_omp) {

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Allgather(&K3_buffer[0], static_cast<int>(n_para_omp*(n_para_mpi/world_size+1)), MPI_DOUBLE,
                  &K3_result[0], static_cast<int>(n_para_omp*(n_para_mpi/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
}

#endif //KELDYSH_MFRG_MPI_SETUP_H
