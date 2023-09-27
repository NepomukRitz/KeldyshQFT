#include <iostream>          // text input/output
#include <sys/stat.h>
#include <bits/stdc++.h>
#include "parameters/master_parameters.hpp"
#include "symmetries/Keldysh_symmetries.hpp"
#include <omp.h>
#include "utilities/mpi_setup.hpp"
#include "mfRG_flow/flow.hpp"
#include "tests/test_perturbation_theory.hpp"
//#include "tests/test_interpolation.hpp"
#include "tests/reproduce_benchmark_data.hpp"
#include "utilities/util.hpp"
#include "utilities/hdf5_routines.hpp"
#include "tests/integrand_tests/saveIntegrand.hpp"
#include "tests/test_symmetries.hpp"
#include "perturbation_theory_and_parquet/perturbation_theory.hpp"
#include "perturbation_theory_and_parquet/parquet_solver.hpp"
#include "tests/test_ODE.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif


auto main(int argc, char * argv[]) -> int {
#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }
#endif
    double t_start = utils::get_time();

    /// Parse and check command line arguments:
    utils::print("number of args: ", argc-1, ", expected: 3 \n");
    const int n_loops = atoi(argv[1]);
    const int n_nodes = atoi(argv[2]);
    //const double U_in = atof(argv[3]);
    const double T_in = atof(argv[3]);
    //const double Gamma_in = atof(argv[5]);
    //const double Vg_in = atof(argv[6]);



    ///fRG runs:
    fRG_config config;
    config.nODE_ = 81;
    config.epsODE_abs_ = 1e-8;
    config.epsODE_rel_ = 1e-6;
    config.nloops = n_loops;
    config.U = 1.0;
    config.T = (ZERO_T ? 0.0 : T_in);
    config.Gamma = ((REG==5) ? 4.0/M_PI : 0.2);
    config.epsilon = (PARTICLE_HOLE_SYMMETRY ? 0.0 : 0.5) - config.U * 0.5;
    config.save_intermediateResults = false;
    config.number_of_nodes = n_nodes;

    utils::check_input(config);
    utils::print_job_info(config);
    std::string filename = utils::generate_filename(config);

    /// Job and Data directory
    std::string job = "T=" + std::to_string(config.T);
    job += "_U=" + std::to_string(config.U);
#ifndef PARTICLE_HOLE_SYMM
    job += "_eVg=" + std::to_string(config.epsilon + config.U*0.5);
#endif
#if SBE_DECOMPOSITION
    job += "_SBE" ;
#endif
    data_dir = utils::generate_data_directory(job);

    if (n_loops > 0){ /// fRG runs:
        n_loop_flow(data_dir+filename, config);
        //get_integrand_dGamma_1Loop<state_datatype>(data_dir, 1, 0);
        //test_PT_state<state_datatype>(data_dir+"sopt.h5", 1.8, false);
    }
    if (n_loops == 0){ /// parquet runs:
        // {5./M_PI*0.5};

        //const std::vector<double> myU_NRG {0.05, 0.25, 0.5, 0.75, 1., 1.25, 1.50, 1.75, 2.0, 2.25, 2.5};
        const std::vector<double> myU_NRG {0.5*M_PI};
        run_parquet(config, myU_NRG, 1, true);
        //run_parquet(config, myU_NRG, 2, true);
        //run_parquet(config, myU_NRG, 3, true);
    }
    if (n_loops == -2){ /// plain and simple second order perturbation theory
        //const std::vector<double> U_over_Delta_list {0.1, 0.05*M_PI, 0.1*M_PI, 0.5, 0.2*M_PI, 0.3*M_PI, 1.0,
        //                                             0.4*M_PI, 1.5, 0.5*M_PI, 0.6*M_PI, 2.0, 0.7*M_PI, 0.75*M_PI, 2.5,
        //                                             0.8*M_PI, 0.9*M_PI, 3.0, 1.0*M_PI, 1.1*M_PI, 3.5, 1.2*M_PI,
        //                                             1.25*M_PI, 4.0, 1.3*M_PI, 1.4*M_PI, 4.5, 1.5*M_PI, 5.0};
        //const std::vector<double> U_over_Delta_list {0.1, 0.2, 0.3, 0.4, 0.5};
        const std::vector<double> U_over_Delta_list{
                0.1     ,  0.109854,  0.120679,  0.132571,  0.145635,  0.15708 ,
                0.159986,  0.175751,  0.19307 ,  0.212095,  0.232995,  0.255955,
                0.281177,  0.308884,  0.314159,  0.339322,  0.372759,  0.409492,
                0.449843,  0.494171,  0.5     ,  0.542868,  0.596362,  0.628319,
                0.655129,  0.719686,  0.790604,  0.868511,  0.942478,  0.954095,
                1.0     ,  1.048113,  1.151395,  1.256637,  1.264855,  1.389495,
                1.5     ,  1.526418,  1.570796,  1.676833,  1.84207 ,  1.884956,
                2.0     ,  2.02359 ,  2.199115,  2.222996,  2.356194,  2.442053,
                2.5     ,  2.513274,  2.682696,  2.827433,  2.947052,  3.0     ,
                3.141593,  3.237458,  3.455752,  3.5     ,  3.55648 ,  3.769911,
                3.90694 ,  3.926991,  4.0     ,  4.08407 ,  4.291934,  4.39823 ,
                4.5     ,  4.712389,  4.714866,  5.0     ,  5.179475,  5.689866,
                6.250552,  6.866488,  7.54312 ,  8.286428,  9.102982, 10.0
        };
        for (double U_over_Delta: U_over_Delta_list) {
            const double Lambda = 2. / U_over_Delta * config.U - config.Gamma;
            State<state_datatype> state (Lambda, config);   // create final and initial state
            state.initialize();
            sopt_state(state);
            const std::string PT2_filename = data_dir + "PT2_U_over_Delta=" + std::to_string(U_over_Delta) +
                    "_T=" + std::to_string(config.T) + "_eVg=" + std::to_string(config.epsilon + config.U*0.5) +
                    "_n1=" + std::to_string(nBOS) + "_n2=" + std::to_string(nBOS2) + "_n3=" + std::to_string(nBOS3) + ".h5";
            write_state_to_hdf(PT2_filename, 0, 1, state);
        }
    }

    //reproduce_benchmark_data();

    utils::print("CPU hours for whole run: \t ");
    utils::get_cpu_hours(t_start);
    utils::hello_world();

#ifdef USE_MPI
    if (MPI_FLAG) {
        MPI_Finalize();
    }
#endif
    return 0;
}
