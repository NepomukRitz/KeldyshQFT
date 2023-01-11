#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "utilities/util.hpp"
#include "utilities/hdf5_routines.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif

auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif
    utils::print(" ---  Post-processing  ---");

    std::string filename = "../Data_KFT=0.100000_U=1.000000_REG=2/parquetInit4_U_over_Delta=0.100000_T=0.100000_eVg=0.000000_n1=201_n2=51_n3=21.h5"; // hdf5 file with states to be postprocessed

    /// functions that do post-processing on several Lambda layers
    compute_Phi_tilde(filename);
    sum_rule_K1tK(filename);
    check_Kramers_Kronig(filename);
    compare_flow_with_FDTs(filename, true);

    /// functions that do postprocessing on a single state:
    const int Lambda_it = 0;    // pick Lambda layer
    State<state_datatype> state = read_state_from_hdf(filename, Lambda_it);
    save_slices_through_fullvertex(filename, 0, 0);


#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}