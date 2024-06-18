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


#include "../utilities/util.hpp"
#include "correlation_functions/state.hpp"
#include "gsl/gsl_interp.h"
#include "build_NRG_state.hpp"

#ifdef USE_MPI
#include <mpi.h>
#endif


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

    /// Code goes here:
    double lambda = 2.0 / U_over_Delta - config.Gamma;

    std::string NRG_DATAPATH = "/Users/nepomuk-work/PhD/NRG_consistency/data/";
    std::string NRG_FILENAME = NRG_DATAPATH + "siam_u1.0.h5";   // TODO: Wrong file for now.

    // TODO: Assertions that the metadata in the NRG file agree with the physical parameters set here.

    // new state to hold NRG data with Hartree value initialized to config.U / 2
    // and vertex initialized to -config.U / 2:
    State<comp,false> NRG_state = State<comp,false>(lambda, config, true);
    NRG_state.vertex.irred().initialize_NRG_input(lambda, config);

    build_NRG_Sigma(NRG_state);
    build_NRG_K1(NRG_state);
    build_NRG_K2_and_K2p(NRG_state);
    build_NRG_rest_term(NRG_state);

    /// Read in data from NRG file
    H5::H5File NRG_file(NRG_FILENAME, H5F_ACC_RDONLY);
    H5::DataSet NRG_dataset = NRG_file.openDataSet("KF/ph/K1/a/up_down/real");
    H5::DataSpace NRG_dataspace = NRG_dataset.getSpace();

    int rank = NRG_dataspace.getSimpleExtentNdims();
    assert(rank==7);
    //utils::print(rank, true);

    hsize_t dims[7];
    NRG_dataspace.getSimpleExtentDims(dims, NULL);

    //utils::print(dims[0], true);
    //utils::print(dims[1], true);
    //utils::print(dims[2], true);
    //utils::print(dims[3], true);
    //utils::print(dims[4], true);
    //utils::print(dims[5], true);
    //utils::print(dims[6], true);

    //NRG_dataset.read();


    std::array<std::size_t, 7> length = {dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], dims[6]};

    multidimensional::multiarray<double,7> K1_a_updown_real_data = multidimensional::multiarray<double,7>(length);
    NRG_dataset.read(K1_a_updown_real_data.data(), H5::PredType::NATIVE_DOUBLE);

    utils::check_input(config);
    utils::hello_world();
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}