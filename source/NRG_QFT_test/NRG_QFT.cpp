/// Main file to perform checks of diagrammatic equations on multipoint NRG data.

/**
* Want to:
 *  - set all parameters as required
 *  - open NRG data from external file
 *  - read in its contents into a self-energy and a vertex (i.e. a state)
 *      - need extra frequency grid?
 *      - need to split up into asymptotic classes?
 *      - Idea: Write R + \sum_r K^{3, r} into the irreducible vertex class.
 *  - plug the vertex into the SDE and generate a new self-energy
 *      - which version to use (I would start with v1)?
 *  - compare the result with the self-energy from NRG
 *      - Don't necessarily need to do this here
*/

#include "../data_structures.hpp"
#include "../utilities/util.hpp"
#include "../utilities/hdf5_routines.hpp"
#include "correlation_functions/state.hpp"
#include "gsl/gsl_interp.h"
#include "read_NRG_data.hpp"
#include "build_NRG_state.hpp"
#include "frequencies_for_NRG.hpp"
#include "identities.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif



State<comp, false> read_or_build_NRG_state(const double& lambda, const fRG_config& config,
                                           const std::string& NRG_FILENAME, const std::string& NRG_Cpp_FILENAME){
    if (std::filesystem::exists(NRG_Cpp_FILENAME)) {
        utils::print("Reading in existing NRG-state ... ");
        State<comp, false> NRG_state = read_state_from_hdf(NRG_Cpp_FILENAME, 0);
        utils::print_add("done.", true);
        return NRG_state;
    }
    else {
        // new state to hold NRG data with Hartree value initialized to config.U / 2
        // and vertex initialized to -config.U / 2:
        State<comp,false> NRG_state = State<comp,false>(lambda, config, true);
        NRG_state.vertex.irred().initialize_NRG_input(lambda, config);

        build_NRG_Sigma(NRG_state, NRG_FILENAME);
        build_NRG_K1(NRG_state, NRG_FILENAME);
        build_NRG_K2_and_K2p(NRG_state, NRG_FILENAME);
        build_NRG_core_as_K3t(NRG_state, NRG_FILENAME);

        write_state_to_hdf(NRG_Cpp_FILENAME, 0, 1, NRG_state);
        return NRG_state;
    }
}


auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif
    /// Parse command line arguments
    const double T_in = atof(argv[1]);              // Temperature in units of U
    const double U_over_Delta = atof(argv[2]);

    /// Parameter assertions
    static_assert(DEBUG_SYMMETRIES == 1);
    static_assert(KELDYSH == true);
    static_assert(ZERO_T == false);
    static_assert(CONTOUR_BASIS == 0);
    static_assert(PARTICLE_HOLE_SYMMETRY == true);
    static_assert(MAX_DIAG_CLASS == 3);
    static_assert(SBE_DECOMPOSITION == 0);
    static_assert(REG == 2);
    static_assert(VECTORIZED_INTEGRATION == 0);

    /// Set up config struct
    fRG_config config;
    config.U = 1.0;
    config.T = T_in;
    config.Gamma = 0.2;
    config.epsilon = - config.U * 0.5;
    config.number_of_nodes = 1;


    const double lambda = 2.0 / U_over_Delta - config.Gamma;

    //std::string NRG_DATAPATH     = "/Users/nepomuk-work/PhD/NRG_consistency/data/";              // for MacBook
    std::string NRG_DATAPATH     = "/dss/dssfs02/pn34vu/pn34vu-dss-0001/ra49hif/mfrg/data/";     // for KCS
    std::string NRG_FILENAME     = NRG_DATAPATH + "siam_u0.5.h5";
    std::string NRG_Cpp_FILENAME = NRG_DATAPATH + "siam_u0.5_C++.h5";

    utils::check_input(config);
    check_NRG_input(NRG_FILENAME, U_over_Delta, T_in);


    /// build required frequency grids to give to MuNRG
    /*
    WantedFrequencyValues freqs = collectWantedFrequencyValues(lambda, config);

    freqs.W_t  = FrequencyProcessorForNRG(freqs.W_t, U_over_Delta).process_frequencies();
    freqs.V_t  = FrequencyProcessorForNRG(freqs.V_t, U_over_Delta).process_frequencies();
    freqs.Vp_t = FrequencyProcessorForNRG(freqs.Vp_t, U_over_Delta).process_frequencies();

    saveWantedFrequenciesToHDF(NRG_DATAPATH + "frequencies.h5", freqs);
    */


    const State<comp, false> NRG_state = read_or_build_NRG_state(lambda, config, NRG_FILENAME, NRG_Cpp_FILENAME);

    const std::string IDENTITIES_FILENAME = NRG_DATAPATH + "siam_u0.5_identities.h5";

    State<comp,false> selfenergy_from_SDE_v3 = evaluate_SDE_from_K1_plus_K2(NRG_state);
    write_state_to_hdf(IDENTITIES_FILENAME, 0, 4, selfenergy_from_SDE_v3);

    State<comp,false> selfenergy_from_SDE_v2 = evaluate_SDE_from_Gamma(NRG_state);
    add_state_to_hdf  (IDENTITIES_FILENAME, 1, selfenergy_from_SDE_v2);

    State<comp,false> K1_from_BSE = evaluate_BSE_for_K1(NRG_state);
    add_state_to_hdf  (IDENTITIES_FILENAME, 2, K1_from_BSE);

    State<comp,false> K1_plus_K2_from_BSE    = evaluate_BSE_for_K1_plus_K2(NRG_state);
    add_state_to_hdf  (IDENTITIES_FILENAME, 3, K1_plus_K2_from_BSE);

    utils::hello_world();
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}