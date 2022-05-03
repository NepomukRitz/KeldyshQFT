#include "mpi_setup.hpp"
#ifdef USE_MPI
int mpi_world_rank() {
    int world_rank= 0;
    if constexpr (MPI_FLAG) MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    return world_rank;
}

int mpi_world_size() {
    int world_size = 1;
    if constexpr (MPI_FLAG) MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    return world_size;
}

void mpi_collect(vec<comp>& buffer, vec<comp>& result, int n_mpi, int n_omp) {
    int world_size = mpi_world_size();
    MPI_Allgather(&buffer[0], static_cast<int>(2*n_omp*(n_mpi/world_size+1)), MPI_COMPLEX,
                  &result[0], static_cast<int>(2*n_omp*(n_mpi/world_size+1)), MPI_COMPLEX, MPI_COMM_WORLD);
}

void mpi_collect(vec<double>& buffer, vec<double>& result, int n_mpi, int n_omp) {
    int world_size = mpi_world_size();

    MPI_Allgather(&buffer[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_DOUBLE,
                  &result[0], static_cast<int>(n_omp*(n_mpi/world_size+1)), MPI_DOUBLE, MPI_COMM_WORLD);
}
#else
int mpi_world_rank() {
    return 0;
}

int mpi_world_size() {
    return 1;
}
void mpi_collect(vec<comp>& buffer, vec<comp>& result, int n_mpi, int n_omp) {
    result = buffer;
}

void mpi_collect(vec<double>& buffer, vec<double>& result, int n_mpi, int n_omp) {
    result = buffer;
}
#endif