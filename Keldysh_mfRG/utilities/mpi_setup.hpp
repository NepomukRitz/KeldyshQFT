/**
 * Functions for handling MPI parallelization: initializing buffers, distributing results between MPI processes etc.
 */

#ifndef KELDYSH_MFRG_MPI_SETUP_HPP
#define KELDYSH_MFRG_MPI_SETUP_HPP

#include "../data_structures.hpp" // real/complex vector classes
#ifdef USE_MPI
#include <mpi.h>             // basic mpi functionality
#endif

// Get the rank(ID) of the current process
int mpi_world_rank();

// Get the number of processes
int mpi_world_size();

/**
 * initialize buffer into which each MPI process can write
 * @tparam Q      : result data type (comp or double)
 * @param n_mpi   : # of MPI tasks (# of grid points for external variable(s) in which we do MPI parallelization)
 * @param n_omp   : # of OMP tasks (# of grid points for external variable(s) in which we do OMP parallelization)
 * @return vec<Q> : empty buffer vector with length corresponding to total number of tasks that the current process
 *                  will compute
 */
template <typename Q>
vec<Q> mpi_initialize_buffer(int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    vec<Q> buffer (n_omp*(n_mpi/world_size+1));
    return buffer;
}

/**
 * initialize result buffer into which results of all MPI processes are collected
 * @tparam Q      : result data type (comp or double)
 * @param n_mpi   : # of MPI tasks (# of grid points for external variable(s) in which we do MPI parallelization)
 * @param n_omp   : # of OMP tasks (# of grid points for external variable(s) in which we do OMP parallelization)
 * @return vec<Q> : empty result buffer vector with length corresponding to the sum of total numbers of tasks
 *                  for all current processes
 */
template <typename Q>
vec<Q> mpi_initialize_result(int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    vec<Q> result (n_omp * (n_mpi - (n_mpi % world_size) + world_size));
    return result;
}

/**
 * distribute results between all MPI processes, using MPI_Allgather (for data type comp)
 * @param buffer : buffer vector containing the results of the current process
 * @param result : buffer vector in which the results of all processes are written
 * @param n_mpi  : # of MPI tasks (# of grid points for external variable(s) in which we do MPI parallelization)
 * @param n_omp  : # of OMP tasks (# of grid points for external variable(s) in which we do OMP parallelization)
 */
void mpi_collect(vec<comp>& buffer, vec<comp>& result, int n_mpi, int n_omp);

/**
 * distribute results between all MPI processes, using MPI_Allgather (for data type double)
 * @param buffer : buffer vector containing the results of the current process
 * @param result : buffer vector in which the results of all processes are written
 * @param n_mpi  : # of MPI tasks (# of grid points for external variable(s) in which we do MPI parallelization)
 * @param n_omp  : # of OMP tasks (# of grid points for external variable(s) in which we do OMP parallelization)
 */
void mpi_collect(vec<double>& buffer, vec<double>& result, int n_mpi, int n_omp);

/**
 * bring results collected from all MPI processes via mpi_collect(...) into the correct order to store them into
 * the vertex objects
 * @tparam Q
 * @param result
 * @param n_mpi
 * @param n_omp
 * @return
 */
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
#endif //KELDYSH_MFRG_MPI_SETUP_HPP
