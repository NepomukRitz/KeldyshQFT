
#include <bits/stdc++.h>
#include <iostream>          // text input/output
#include <sys/stat.h>
#include "../../parameters/master_parameters.hpp"
#include "saveIntegrand.hpp"

auto main(int argc, char * argv[]) -> int {
    std::string job = "U=" + std::to_string(glb_U);
    data_dir = "../Data_KF" + job + "/";
    std::string dir_str, filename;
    char channel;
    int it_Lambda, k_class_int, rkStep, i0, i2, i_in;
    double w, v, vp;
    std::cout << "----  Getting integrand  ----" << std::endl;
    std::cout << "number of args: " << argc-1 << ", expected: 11" << std::endl;
    /// Prompt user for input:
    // dir_str: directory in which the intermediate results lie
    // it_Lambda: iteration of ODE solver
    // k_class_int: integer for the K_class 0->k1, 1->k2, 2->k2', 3->k3
    // channel: char for channel
    // i0: external Keldysh indices ranging in [0,...,15]
    // i2: internal Keldysh indices ranging in [0,..., 9] (--> directly corresponding to non-zero components of the BubbleObject)
    // w, v, vp: frequencies in the natural parametrization of channel
    // i_in: internal index
    //std::cin >> dir_str >> it_Lambda >> k_class_int >> channel >> i0 >> i2 >> w >> v >> vp >> i_in;
    /// Parse input:
    dir_str = argv[1];
    filename = argv[2];
    it_Lambda = atoi(argv[2]);
    rkStep = atoi(argv[3]);
    k_class_int = atoi(argv[4]);
    channel = *(argv[5]);
    i0 = atoi(argv[6]);
    i2 = atoi(argv[7]);
    int spin = atoi(argv[8]);
    w = atof(argv[9]);
    v = atof(argv[10]);
    vp = atof(argv[11]);
    i_in = atoi(argv[12]);
    K_class k_class = static_cast<K_class>(k_class_int);

    /// print input arguments:
    std::cout << "Check the input arguments: " << std::endl;
    std::cout << "dir_str: " << dir_str << ", it_Lambda: " << it_Lambda << ", k_class_int: " << k_class_int
    << ", channel: " << channel << ", i0: " << i0 << ", i2: " << i2
    << ", w: " << w << ", v: " << v << ", vp: " << vp << ", i_in: " << i_in << std::endl;

    dir_str = dir_str + "intermediateResults/";
    const std::string file_Psi = dir_str + filename;


    std::string dir_integrand_str = "integrands/";
    makedir(data_dir + dir_integrand_str);
    const std::string filename_prefix = dir_integrand_str + "dGamma1Loop_iLambda"+std::to_string(it_Lambda)+"_RKstep"+std::to_string(rkStep);


    saveIntegrand::dGamma_1Loop_fromFlow<state_datatype>(filename_prefix, file_Psi, it_Lambda, k_class, channel, i0, i2, w, v, vp, i_in);
    std::cout << "Integrand for dGamma1_loop successfully created." << std::endl;

    return 0;

}