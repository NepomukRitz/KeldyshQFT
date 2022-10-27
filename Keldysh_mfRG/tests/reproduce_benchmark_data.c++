#include "reproduce_benchmark_data.hpp"

void reproduce_benchmark_data() {

    std::cout << " --- Reproducing benchmark data --- " << std::endl;

    const double Lambda = 2.;
    fRG_config config;
    config.nODE_ = 1;
    config.epsODE_abs_ = 1e-8;
    config.epsODE_rel_ = 1e-5;
    config.nloops = 4;
    config.U = 1.0;
    config.T = (ZERO_T ? 0.0 : 0.1);
    config.Gamma = 0.2;
    config.epsilon = (PARTICLE_HOLE_SYMMETRY ? 0.0 : 0.5) - config.U * 0.5;
    config.save_intermediateResults = true;
    config.number_of_nodes = 1;
    utils::check_input(config);
    utils::print_job_info(config);


    State<state_datatype> state_ini(Lambda, config);   // create final and initial state
    state_ini.initialize();     // initialize state with bare vertex and Hartree term in selfenergy
    sopt_state(state_ini);

    const std::string bm_data_dir = "../BmDat_";
    data_dir = bm_data_dir + (KELDYSH_FORMALISM ? "KF" : "MF") + "_" + (PARTICLE_HOLE_SYMMETRY ? "PHS" : "nonPHS") + "_" + (ZERO_T ? "T0" : "T1e-2") + "_" + (SBE_DECOMPOSITION ? "SBE" : "ASY") + "_REG" + std::to_string(REG) + "/";
    utils::makedir(data_dir);

    std::string parquet_filename_withK3 = data_dir + "parquetInit4_final_n1=" + std::to_string(nBOS) + (MAX_DIAG_CLASS > 1 or true? "_n2=" + std::to_string(nBOS2) + (MAX_DIAG_CLASS > 2 or true? "_n3=" + std::to_string(nBOS3) : "") : "") + ".h5";
    parquet_solver(parquet_filename_withK3, state_ini, Lambda, 2, 1e-4, 2  );


    State<state_datatype> dPsi_dLambda(Lambda, config);
    rhs_n_loop_flow_t<state_datatype> rhs(config);
    rhs(state_ini, dPsi_dLambda, Lambda);


}

