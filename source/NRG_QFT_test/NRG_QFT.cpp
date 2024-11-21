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

#include <optional>
#include <filesystem>

#ifdef USE_MPI
#include <mpi.h>
#endif



State<comp, false> read_or_build_NRG_state(const double& lambda, const fRG_config& config,
                                           const std::string& MuNRG_FILENAME, const std::string& NRG_Cpp_FILENAME,
                                           const std::optional<std::string>& NRG_SELFENERGY_FILENAME = std::nullopt){
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

        if (NRG_SELFENERGY_FILENAME.has_value()) {
            build_NRG_Sigma(NRG_state, NRG_SELFENERGY_FILENAME.value());
        }
        else {
            build_NRG_Sigma_from_MuNRG(NRG_state, MuNRG_FILENAME);
        }

        build_NRG_K1(NRG_state, MuNRG_FILENAME);
        build_NRG_K2_and_K2p(NRG_state, MuNRG_FILENAME);
        build_NRG_core_as_K3t(NRG_state, MuNRG_FILENAME);

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
    const double u = atof(argv[2]);                 // value for u = U / (πΔ)
    const double U_over_Delta = u * M_PI;
    const double D_in = atof(argv[3]);              // hybridization band-width. Use >= 10000 for wide-band limit

    /// Parameter assertions
    static_assert(DEBUG_SYMMETRIES == 1);
    static_assert(KELDYSH == true);
    static_assert(ZERO_T == false);
    static_assert(CONTOUR_BASIS == 0);
    static_assert(PARTICLE_HOLE_SYMMETRY == true);
    static_assert(MAX_DIAG_CLASS == 3);
    static_assert(SBE_DECOMPOSITION == 0);
    static_assert(REG == 2);
    static_assert(VECTORIZED_INTEGRATION == 1);

    /// Set up config struct
    fRG_config config;
    config.U = 1.0;
    config.T = T_in;
    config.Gamma = 0.2;
    config.epsilon = - config.U * 0.5;
    config.D = D_in;
    config.number_of_nodes = 1;


    const double lambda = 2.0 / U_over_Delta - config.Gamma;

    std::ostringstream u_str;
    u_str << std::fixed << std::setprecision(1) << u;

    //const std::string NRG_DATAPATH        = "/Users/nepomuk-work/PhD/NRG_consistency/data/weak/";              // for MacBook
    const std::string NRG_DATAPATH        = "/dss/dssfs02/pn34vu/pn34vu-dss-0001/ra49hif/mfrg/data/";     // for KCS
    const std::string MuNRG_FILENAME      = NRG_DATAPATH + "siam_u"+u_str.str()+".h5";
    const std::string NRG_Cpp_FILENAME    = NRG_DATAPATH + "siam_u"+u_str.str()+"_C++.h5";
    const std::string IDENTITIES_FILENAME = NRG_DATAPATH + "siam_u"+u_str.str()+"_identities.h5";

    const std::string NRG_SELFENERGY_FILENAME = NRG_DATAPATH + "SIAM_NRG4fRG_Gamma=1_U=1.5708_T=0.015708_eVg=0_Lambda=2_nz=6_Nkeep=5000_Etrunc=12.h5"; // only weak coupling as of now

    utils::check_input(config);
    check_NRG_input(MuNRG_FILENAME, U_over_Delta, T_in);


    /// build required frequency grids to give to MuNRG
    /*
    WantedFrequencyValues freqs = collectWantedFrequencyValues(lambda, config);

    freqs.W_t  = FrequencyProcessorForNRG(freqs.W_t, U_over_Delta).process_frequencies();
    freqs.V_t  = FrequencyProcessorForNRG(freqs.V_t, U_over_Delta).process_frequencies();
    freqs.Vp_t = FrequencyProcessorForNRG(freqs.Vp_t, U_over_Delta).process_frequencies();

    saveWantedFrequenciesToHDF(NRG_DATAPATH + "frequencies.h5", freqs);
    */


    const State<comp, false> NRG_state = read_or_build_NRG_state(lambda, config,
                                                                 MuNRG_FILENAME,
                                                                 NRG_Cpp_FILENAME,
                                                                 NRG_SELFENERGY_FILENAME);


    IdentityChecker Identities(NRG_state, NRG_DATAPATH + "siam_u"+u_str.str());
    Identities.check_BSE();
    Identities.check_SDE();
    Identities.compute_1D_WardIdentity_wrt_v_RHS();
    Identities.compute_1D_WardIdentity_wrt_w_RHS();
    Identities.compute_2D_WardIdentity();


    utils::hello_world();
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}