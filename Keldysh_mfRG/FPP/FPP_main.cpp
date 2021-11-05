//
// Created by Marcel on 14.04.2021.
//

//#define INTEGRATOR_TYPE 5;

#include <iostream>
#include <iomanip>
#include <bits/stdc++.h>
#include "Coordinates.h"
#include "Momentum-integral-Bubble.h"
#include "Differential_Equation.h"
#include "Ladder-approximation.h"
#include "selfenergy_loop.h"
#include "fRG-T-matrix-approach.h"
#include "1D-integrals.h"
#include "../utilities/util.h"
#include "../ODE_solvers.h"
#include "../grids/flow_grid.h"
#include "../data_structures.h"
//#include "../OldFiles/paid.hpp"
//#include "zeros.h"


int main() {

    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }

    gsl_set_error_handler_off();

    std::cout << "Hello world!\n";
    std::cout << "test \n";

    double t0 = get_time();

    // Coordinate Transformations
    // =======================================

    /*double x;
    double y;
    std::cout << "Type in x.\n";
    cin >> x;
    std::cout << "Type in y.\n";
    cin >> y;
    double angle;
    angle = phi(x, y);
    double angledegree = angle*180/M_PI;
    std::cout << "Phi is " << angledegree << ".\n";
    return 0;*/
    //coordinate_transform();
    //coordinate_transform();
    /*vector<double> v_in;
    v_in = {1.0,1.0,1.0};
    vector<double> v_out;
    v_out = cartesion_to_spherical(v_in);
    std::cout << v_out << ".\n"; */
    /* for (int n = 1; n < 100; ++n) {
        vector<double> v_in = {1.1223498, -2.24342338, 3.322433};
        vector<double> v_out;
        v_out = transform_back_forth(v_in, n);
        double vx = v_out[0];
        double vy = v_out[1];
        double vz = v_out[2];
        std::cout << "N = " << n << ": v = {" << setprecision(8) << vx << "," << setprecision(8) << vy << "," << setprecision(8) << vz << "}.\n";
    } */

    //monte_carlo_test1();

    // Parameters bare Green's function
    // =======================================

    glb_muc = 1.0;
    glb_mud = 0.0;
    glb_mc = 1.0;
    glb_md = 1.0;
    glb_ainv = 1.0;
    glb_prec = 0.;

    // std::cout Bubble
    // =======================================

    /*std::cout << "Give nu, k^2, i\n";
    double nu;
    double ksquared;
    char i;
    cin >> nu;
    cin >> ksquared;
    cin >> i;
    std::cout << "Your bare Green's function is :" << G0(nu, ksquared, i) << "\n";
    std::cout << "Give omega, nu, k, q, cos(theta), channel, i, j\n";
    double omega;
    double k;
    double q;
    double x;
    char channel;
    char j;
    cin >> omega;
    cin >> nu;
    cin >> k;
    cin >> q;
    cin >> x;
    cin >> channel;
    cin >> i;
    cin >> j;
    std::cout << "Your bare bubble is: " << SimpleBubble(omega, nu, q, k, x, channel, i, j) << "\n"; */

    // Calculate Bubble integrals
    // =======================================

    //integrate_bubble_full_monte_carlo(0.0,0.0,0.0,'c','d','p',1000.,50000,0.1);
    //comp output = integrate_bubble_vegas(0.0,0.0,0.0,'c','d','p',1000.,50000,0.1);
    //std::cout << "Result: " << output << "\n";

    //IMPORTANT LINE IS FOLLOWING
    // integral_bubble_w_vpp_list ('c', 'c', 'p', 2.0, 2.0, 1.0, 1000., 50000, 0.1, 40, 40, 5);

    //Monte Carlo Test with parameters
    // =======================================
    /*
    gsl_monte_function F;
    struct my_f_params params = {3.0, 2.0, 1.0};

    F.f = &my_f;
    F.dim = 2;
    F.params = &params;

    double res, err;

    const gsl_rng_type *T;
    gsl_rng *r;

    double xl[2] = { 0., 0.};
    double xu[2] = { 1., 1.};

    size_t calls = 50000;

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
        gsl_monte_plain_integrate (&F, xl, xu, 2, calls, r, s,
                                   &res, &err);

        gsl_monte_plain_free (s);

        std::cout << "result: " << res << ", error: " << err << "\n";
    }
    //*/

    // Test if clause with char
    //====================================
    /*
    char c;
    std::cout << "Type in character: ";
    cin >> c;
    double output2 = test2(c);
    std::cout << "\n Your value is " << output2 << ".\n";
    */

    // 1D INTEGRALS GSL
    // ===================================

    //integral_1D_with_singularity();
    //integral_1D_with_infinite_range();

    /*
    // Exact bubble and numerical bubble
    print_exact_bubble(3.4,-12.3,0.0,'c','d','a');
    print_numerical_bubble(3.4,-12.3,0.0,'c','d','a',0.,0.);
    print_exact_bubble(3.4,-12.3,0.0,'c','d','p');
    print_numerical_bubble(3.4,-12.3,0.0,'c','d','p',0.,0.);
    print_exact_bubble(3.4,-12.3,0.0,'c','d','t');
    print_numerical_bubble(3.4,-12.3,0.0,'c','d','t',0.,0.);
    print_exact_bubble(3.4,-12.3,1.0,'c','d','a');
    print_numerical_bubble(3.4,-12.3,1.0,'c','d','a',0.,0.);
    print_exact_bubble(3.4,-12.3,1.0,'c','d','p');
    print_numerical_bubble(3.4,-12.3,1.0,'c','d','p',0.,0.);
    print_exact_bubble(3.4,-12.3,1.0,'c','d','t');
    print_numerical_bubble(3.4,-12.3,1.0,'c','d','t',0.,0.);
    comp test001 = -1/(2*sqrt(2)*M_PI*sqrt(-1));
    std::cout << test001 << "\n";
    comp test002 = pow(-3.0+glb_i,0.5);
    std::cout << test002 << "\n";
    comp test003 = -1./(2*sqrt(2)*M_PI*pow(-1.0+1e-16*glb_i,0.5));
    std::cout << test003 << "\n";
    print_exact_bubble(0.0,0.0,0.0,'c','d','p');
    print_numerical_bubble(0.1,0.1,0.1,'c','d','p',0.,0.);
    print_numerical_bubble(0.,0.,0.,'c','d','p',0.,0.);


    //integral_bubble_w_vpp_list_exact ('c', 'c', 'p', 10., 10., 10., 200, 200, 11);
    integral_bubble_w_vpp_list_exact ('c', 'c', 'p', 1., 2., 3., 4, 5, 6);
    //integral_bubble_w_vpp_list_exact ('c', 'c', 't', 10., 10., 10., 200, 200, 11);
    //integral_bubble_w_vpp_list_exact ('c', 'c', 'a', 10., 10., 10., 200, 200, 11);

    //integral_bubble_w_vpp_list_2D ('c', 'c', 'p', 10., 10., 10., 200, 200, 11,10e-12,10e-12);
    //integral_bubble_w_vpp_list_2D ('c', 'c', 't', 10., 10., 10., 200, 200, 11,10e-12,10e-12);
    //integral_bubble_w_vpp_list_2D ('c', 'c', 'a', 10., 10., 10., 200, 200, 11,10e-12,10e-12);
    */

    /*double output1, output2;
    comp output3;
    output1 = bubble_integrate_theta (0.1, 0.1, 0.1, 0.1, 'c', 'c', 0, 0.);
    std::cout << "Output 1: " << output1 << "\n";
    output2 = bubble_integrate_kpp (0.1, 0.1, 0.1, 'c', 'c', 0, 0.,0.);
    std::cout << "Output 2: " << output2 << "\n";
    output3 = bubble_k_2d_integrated (0.1, 0.1, 0.1, 'c', 'c', 'p', 0, 0);
    double reoutput3, imoutput3;
    reoutput3 = real(output3);
    imoutput3 = imag(output3);
    std::cout << "Output 3: " << reoutput3 << " + i " << imoutput3 << "\n";
    print_numerical_bubble(0.1, 0.1, 0.1, 'c', 'c', 'p', 0, 0);*/


    //std::cout << "now new" << "\n";

    //print_numerical_bubble(-2.0,2.0,10.,'c','c','p',0,0);
    //print_numerical_bubble(-4.0,-4.0,8.0,'c','d','t',0,0);
    //print_numerical_bubble(0.0,0.0,0.0,'c','d','t',0,0);

    // Loop_Integral

    /*print_numerical_loop (0.1, 0.1, 'c', 0);
    print_numerical_loop (0.01, 0.1, 'c', 0);
    print_numerical_loop (0.001, 0.1, 'c', 0);*/

    //integral_loop_Lambda_vp_list ('c', 0.1, 100, 10, 100, 200, 0);

    // Differential equations
    // =============================

    /*double test_heaviside, x_pos, x_neg, x_zero;
    x_pos = 353.2;
    x_neg = -1.3e-16;
    x_zero = 0.000000;
    test_heaviside = heaviside(x_pos);
    std::cout << test_heaviside << "\n";
    test_heaviside = heaviside(x_neg);
    std::cout << test_heaviside << "\n";
    test_heaviside = heaviside(x_zero);
    std::cout << test_heaviside << "\n";*/


    /*comp test_sfbb_a, test_sfbb_p, test_sfbb_t;
    test_sfbb_a = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 'a');
    std::cout << test_sfbb_a << "\n";
    test_sfbb_p = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 'p');
    std::cout << test_sfbb_p << "\n";
    test_sfbb_t = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 't');
    std::cout << test_sfbb_t << "\n";

    std::cout << "now solving ODE: \n";
    double w = 3.4, q = 1.2, g = -0.5, Lambda_i = 100.0, Lambda_f = 0.001;
    comp outputK1p = K1cdcd_solution(w, q, g, Lambda_i, 'p',Lambda_f, 1e-10, 1e-10, 0.0);
    comp outputK1a = K1cdcd_solution(w, q, g, Lambda_i, 'a',Lambda_f, 1e-10, 1e-10, 0.0);
    comp outputK1t = K1cdcd_solution(w, q, g, Lambda_i, 't',Lambda_f, 1e-10, 1e-10, 0.0);
    std::cout << "K1p = " << outputK1p << "\n";
    std::cout << "K1a = " << outputK1a << "\n";
    std::cout << "K1t = " << outputK1t << "\n";
    w = 0.0, q = 0.0, g = -0.1, Lambda_i = 100.0, Lambda_f = 1e-10;
    outputK1p = K1cdcd_solution(w, q, g, Lambda_i, 'p',Lambda_f, 1e-10, 1e-10, 0.0);
    outputK1a = K1cdcd_solution(w, q, g, Lambda_i, 'a',Lambda_f, 1e-10, 1e-10, 0.0);
    outputK1t = K1cdcd_solution(w, q, g, Lambda_i, 't',Lambda_f, 1e-10, 1e-10, 0.0);
    std::cout << "K1p = " << outputK1p << "\n";
    std::cout << "K1a = " << outputK1a << "\n";
    std::cout << "K1t = " << outputK1t << "\n";
    *//*outputK1p = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1a = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1t = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    std::cout << "K1p = " << outputK1p << "\n";
    std::cout << "K1a = " << outputK1a << "\n";
    std::cout << "K1t = " << outputK1t << "\n";*/
    /*
    std::cout << "now exact with w = 0 = q: \n";
    outputK1p = K1cdcd_solution(0., 0., -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0);
    outputK1a = K1cdcd_solution(0., 0., -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0);
    outputK1t = K1cdcd_solution(0., 0., -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0);
    std::cout << "K1p = " << outputK1p << "\n";
    std::cout << "K1a = " << outputK1a << "\n";
    std::cout << "K1t = " << outputK1t << "\n";
    //K1Lambda (0.0, 0.0, -0.1, 'p', 10.0, 1e-10, 100, 1e-12, 1e-12, 0.0);
    //K1Lambdag(0.0,0.0,-0.1,0.1,'p',10.0,1e-10,100,1e-12,1e-12,0.0);
    *//*outputK1p = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1a = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1t = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    std::cout << "K1p = " << outputK1p << "\n";
    std::cout << "K1a = " << outputK1a << "\n";
    std::cout << "K1t = " << outputK1t << "\n";*/

    //K1Lambda (0.1, 0.1, -0.1, 'p', 10.0, 1e-10, 100, 1e-12, 1e-12, 0.0, 1e-12, 1e-12);

        /*derivativefsqrt();
        ODE_solver_test();
        double outputre, outputim;
        comp outputcomp;
        outputre = dLsfebb ( 0.1, 1.0, 0.1, 'c', 'c', 'p', 0, 1e-10);
        outputim = dLsfebb ( 0.1, 1.0, 0.1, 'c', 'c', 'p', 1, 1e-10);
        outputcomp = outputre + glb_i * outputim;
        std::cout << "Derivative is " << outputcomp << "\n";
        ODE_solver_K1p(3.0, 3.0, 0.1, 10.0, -0.1, 10);*/

    // Sharp frequency regulator
    // =============================

    /*comp output13 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'a');
    comp output14 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'p');
    comp output15 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'t');
    std::cout << "Sharp regulator: " << output13 << ".\n";
    std::cout << "Sharp regulator: " << output14 << ".\n";
    std::cout << "Sharp regulator: " << output15 << ".\n";*/


    // DOUBLE PENDULUM
    // ===============================

    /*std::cout << "now double pendulum: " << "\n";
    solve_double_pendulum(9.81,1.0,2.0,1.0,1.5,0.0,10.0,0.0,30.0,1.0,-2.3);*/

    // ZEROS
    // ==============================
    /*
    double testroot1 = find_root_divergence (sigmoidal, -1.0, 1., 1000, 10e-10);
    std::cout << "x = " << testroot1 << "\n";
    comp exactvalue = sigmoidal(testroot1);
    std::cout << "sigmoid(x) = " << exactvalue << "\n";

    double testroot2 = fine_root_newton (fzeros1, 3, 10e-10, 1000, 10e-10);
    std::cout << "x = " << testroot2 << "\n";
    double exactvalue2 = fzeros1(testroot2);
    std::cout << "f1(x) = " << exactvalue2 << "\n";

    double testroot3 = fine_root_newton (fzeros2, -1.0, 1., 1000, 10e-10);
    std::cout << "x = " << testroot3 << "\n";
    double xx0 = pow(M_PI,3)/log(2)/(10*M_EULER);
    std::cout << "x_exact " << xx0 << "\n";
    double exactvalue3 = fzeros2(testroot3);
    std::cout << "f2(x) = " << exactvalue3 << "\n";
    */

    // LADDER APPROXIMATION
    // ==============================

    glb_muc = 0.0;
    glb_mud = 0.0;
    glb_ainv = 1.0;

    // test sort algorithm

    std::cout << "test sort algorithm:\n";
    double Lambda_i = 1e4;
    double Lambda_f = 1e-10;
    double Lambdaif = 1000.;
    double Lambda_one = 1.0;
    double mylist[] = {1.0, Lambda_i*Lambda_f};
    std::vector<double> myvector (mylist, mylist+2);               // 32 71 12 45 26 80 53 33
    std::cout << "finished\n";
    // using default comparison (operator <):
    std::sort (myvector.begin(), myvector.begin()+2);
    std::cout << "sort finished \n";
    std::cout << "0th: " << myvector[0] << "\n";
    std::cout << "1st: " << myvector[1] << "\n";
    /*for (int i = 0; i < 2; ++i){
        std::cout << myvector[i] << ", ";
    }*/


    /*
    comp bubble1 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 10, 1e-10, 0);
    double bubble1e = exactzerobubble(10, 1e-10);
    //comp output2t = perform_vacuum_integral (10, 1e-10);
    comp bubble2 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 1e4, 1e-10, 0);
    double bubble2e = exactzerobubble(1e4, 1e-10);
    //comp output3t = perform_vacuum_integral (1, 1, 1e4, 1e-10);
    comp bubble3 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 10, 1e-5,0);
    double bubble3e = exactzerobubble(10, 1e-5);
    //comp output4t = perform_vacuum_integral (1, 1, 10, 1e-5);
    comp bubble4 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 1e4, 1e-5,0);
    double bubble4e = exactzerobubble(1e4, 1e-5);
    //comp output5t = perform_vacuum_integral (1, 1, 1e4, 1e-5); */
    /*
    std::cout << "Bubble integral = " << bubble1 << ", exact = " << bubble1e << "\n";
    //double test1 = 1e4;
    //double test2 = sqrt(test1);
    //std::cout << test1 << " = " << test2 << "\n";
    std::cout << "Bubble integral = " << bubble2 << ", exact = " << bubble2e << "\n";
    std::cout << "Bubble integral = " << bubble3 << ", exact = " << bubble3e << "\n";
    std::cout << "Bubble integral = " << bubble4 << ", exact = " << bubble4e << "\n";
    std::cout << "test now \n";
    glb_muc = 0.0;
    glb_mud = 0.0;
    glb_ainv = 1.0;

    comp ladder1, ladder2, ladder3, ladder4, ladder5;
    double laddere;
    std::cout << "k integral exact\n";
    ladder1 = ladder(0.0,0.0,'p',10,1e-10,1,0);
    ladder2 = ladder(0.0,0.0,'p',1e4,1e-10,1,0);
    ladder3 = ladder(0.0,0.0,'p',10,1e-5,1,0);
    ladder4 = ladder(0.0,0.0,'p',1e2,1e-5,1,0);
    ladder5 = ladder(0.0,0.0,'p',1e4,1e-5,1,0);
    laddere = -4*M_PI;
    std::cout << "ladder = " << ladder1 << "\n";
    std::cout << "ladder = " << ladder2 << "\n";
    std::cout << "ladder = " << ladder3 << "\n";
    std::cout << "ladder = " << ladder4 << "\n";
    std::cout << "ladder = " << ladder5 << "\n";
    std::cout << "exact = " << laddere << "\n";



    std::cout << "k integral numerical\n";
    ladder1 = ladder(0.0,0.0,'p',10,1e-10,1,1);
    ladder2 = ladder(0.0,0.0,'p',1e4,1e-10,1,1);
    ladder3 = ladder(0.0,0.0,'p',10,1e-5,1,1);
    ladder4 = ladder(0.0,0.0,'p',1e2,1e-5,1,1);
    ladder5 = ladder(0.0,0.0,'p',1e4,1e-5,1,1);
    laddere = -4*M_PI;
    std::cout << "ladder = " << ladder1 << "\n";
    std::cout << "ladder = " << ladder2 << "\n";
    std::cout << "ladder = " << ladder3 << "\n";
    std::cout << "ladder = " << ladder4 << "\n";
    std::cout << "ladder = " << ladder5 << "\n";
    std::cout << "exact = " << laddere << "\n";

    glb_muc = 1.0;
    glb_ainv = 1.0;
    glb_mud = -1.24;
    comp ladder01 = 1./ladder(0.0,0.0,'p',1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    glb_mud = -1.23;
    ladder01 = 1./ladder(0.0,0.0,'p',1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    glb_mud = -1.2325;
    ladder01 = 1./ladder(0.0,0.0,'p',1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    double root_test0 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,-2.,1e-10,100,1e-16);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test0 << "\n";


    glb_ainv = 2.0;
    double root_test1 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,-2.,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test1 << "\n";
    glb_ainv = 1.5;
    double root_test2 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test1,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test2 << "\n";
    glb_ainv = 1.0;
    double root_test3 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test2,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test3 << "\n";
    glb_ainv = 0.5;
    double root_test4 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test3,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test4 << "\n";
    glb_ainv = 0.0;
    double root_test5 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test4,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test5 << "\n";
    glb_ainv = -0.5;
    double root_test6 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test5,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test6 << "\n";
    glb_ainv = -1.0;
    double root_test7 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test6,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test7 << "\n";
    glb_ainv = -2.0;
    double root_test8 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,root_test7,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test8 << "\n";
    glb_ainv = 0.282843;
    double root_test9 = find_root_newton(0.0,0.0,'p',1e4,1e-10,1,0,-0.342722,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test9 << "\n";

    std::cout << "now solver with exact momentum integral:\n";
    ladder_list ('p',1e4, 1e-10, 1, 0, -2.*sqrt(2), 2.*sqrt(2), -10., 1e-10, 100, 1e-16, 301);
    std::cout << "now solver with numerical momentum integral:\n";
    //ladder_list ('p',1e4,1e-10,1,1,-2*sqrt(2),2*sqrt(2),-10.,1e-10,100,1e-16,101);
    */

    // comp testfmu = test_f_mu(0.0,0.0,1e3,1e-4,1,1.0);
    // std::cout << "f_mu = " << testfmu << "\n";

    // 1-LOOP FRG
    // ==================================
    /*
    glb_ainv = 1.0;
    glb_muc = 1.0;
    glb_mud = -5.;
    std::cout << "bare bubble with numerical k-integral:\n";
    comp testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','t',0);
    std::cout << "t: " << testrhs << "\n";
    testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','p',0);
    std::cout << "p: " << testrhs << "\n";
    testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','a',0);
    std::cout << "a: " << testrhs << "\n";

    std::cout << "bare bubble with exact k-integral:\n";
    testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','t',1);
    std::cout << "t: " << testrhs << "\n";
    testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','p',1);
    std::cout << "p: " << testrhs << "\n";
    testrhs = sharp_frequency_bare_bubble(0.,10.,0.1,'c','d','a',1);
    std::cout << "a: " << testrhs << "\n";
     */
    /*
    comp Gamma_fRG_exact = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-10, 1, 0);
    comp Gamma_fRG_num = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-10, 1, 1);
    comp Gamma_ladder_exact = ladder(0.1,0.1,1e4,1e-10,1,0);
    comp Gamma_ladder_num = ladder(0.1,0.1,1e4,1e-10,1,1);

    std::cout << "Gamma_fRG exact k-integral: " << Gamma_fRG_exact << "\n";
    std::cout << "Gamma_fRG numerical k-integral: " << Gamma_fRG_num << "\n";
    std::cout << "Gamma_ladder exact k-integral: " << Gamma_ladder_exact << "\n";
    std::cout << "Gamma_ladder numerical k-integral: " << Gamma_ladder_num << "\n";
    */
    /*
    comp yy_fin;
    const comp yy_ini = 0.0; // (11.3834,-1.14485e-15) + gint(1e4,1e-10,1);
    const double L_ini = 1e4;
    const double L_fin = 1e-10;

    std::cout << "solve 1lfRG" << "\n";
    ODE_solver_RK4(yy_fin, L_fin, yy_ini, L_ini,
                   rhs_test2, log_substitution, log_resubstitution, nODE);
    */

    // comp ladder_compare = ladder(0.1,0.1,1e4,1e-10,1);
    //std::cout << "yfin = " << Gamma_fRG << ", ladder = " << ladder_compare << "\n";
    //std::cout << "nODE = " << nODE << "\n";

    glb_muc = 1.0;
    glb_ainv = 2.0;
    // double mud_before = glb_mud;
    //double mufRG1 = find_root_fRG (0.0,0.0,1e4,1e-10,1,-5.0,0.5,0.5,100,1e-16);
    //double muladder1 = find_root_ladder (0.0,0.0,1e4,1e-10,1,-5.0,0.5,1e-10,100,1e-16);
    //std::cout << "a^{-1} = " << glb_ainv << ": mu_fRG = " << mufRG1 << ", mu_ladder = " << muladder1 << "\n";
    //std::cout << "mud_i = " << mud_before << ", mud_f = " << glb_mud << "\n";

    //fRG_p_list (1e4, 1e-10, 1,0, -2*sqrt(2), 2*sqrt(2), -10.0, 0.5, 100, 1e-16, 101);
    //ladder_p_list (1e4, 1e-10, 1,0,-2*sqrt(2), 2*sqrt(2), -10.0, 1e-10, 100, 1e-16, 301);
    // fRG_p_list (1e4, 1e-10, 1,1, -2*sqrt(2), 2*sqrt(2), -10.0, 0.5, 100, 1e-16, 11);

    glb_muc = 1.0;
    glb_mud = 0.0;
    glb_ainv = 0.0;

    /* finding divergences in bare bubble
    double gint_test = gint(1e4,1e-10,1);
    std::cout << "g test = " << gint_test << "\n";
    comp limit_test = exact_bare_bubble (-2.0, -1., 0., 'c', 'd', 'p');
    std::cout << "limit test = " << limit_test << "\n";
    std::cout << "abs limit test = " << std::abs(limit_test) << "\n";
    if (isnan(std::abs(exact_bare_bubble (-2.0, -1., 0., 'c', 'd', 'p')))){
        std::cout << "too large \n";
    }
    limit_test = exact_bare_bubble (-2.0, -1.+1e-16, 0., 'c', 'd', 'p');
    std::cout << "limit test = " << limit_test << "\n";
    limit_test = exact_bare_bubble (-2.0, -1.-1e-10, 0., 'c', 'd', 'p');
    std::cout << "limit test = " << limit_test << "\n";
    comp test_comp_sqrt = pow(-1./(2.*(1.+glb_i*(-2.))),-0.5);
    std::cout << "sqrt1 = " << test_comp_sqrt << "\n";
    test_comp_sqrt = pow(-1./(2.*(glb_i*(-0.))),-0.5);
    std::cout << "sqrt2 = " << test_comp_sqrt << "\n";
    comp bubble_test = perform_Pi0_vpp_integral (-2.0, 0., 'c', 'd', 'p', 1e4, 1e-10);
    std::cout << "bubble test = " << bubble_test << "\n";
    comp testm50 = ladder(-2.0, 0., 1e4, 1e-10, 1);

    std::cout << "test = " << testm50 << "\n"; */


    double w, q;
    comp testfRGwq;
    vec<double> ws(6), qs(6);
    vec<comp> testfRGwqs(36);

    /*
    testfRGwq = fRG_solve_nsc(-10., 10., 1e4, 1e-10, 1);
    std::cout << "Gamma_fRG = " << testfRGwq << "\n";


    for (int wi = 0; wi < 6; ++wi){
        w = -10. + wi*std::abs(20.)/5.;
        ws[wi] = w;
        for (int qi = 0; qi < 6; ++qi){
            q = qi*std::abs(10.)/5.;
            qs[qi] = q;
            testfRGwq = fRG_solve_nsc(w, q, 1e4, 1e-10, 1);
            testfRGwqs[wi*6+qi] = testfRGwq;
            //std::cout << "w = " << w << ", q = " << q << ": result = " << testladderwq << "\n";
        }
    }
    std::cout << "now list test fRG: \n";
    for (int wi = 0; wi < 6; ++wi){
        for (int qi = 0; qi <6; ++qi){
            std::cout << "w = " << ws[wi] << ", q = " << qs[qi] << ", Gamma = " <<  testfRGwqs[wi*6+qi] << "\n";
        }
    }
    */

    double t_par = get_time();
    /*
    std::cout << "now list fRG: \n";
    fRG_list_wq(10.,10.,'p',1e4,1e-10,1,0,201,101);
    fRG_list_wq(10.,10.,'a',1e4,1e-10,1,0,201,101);
    get_time(t_par);
    std::cout << "now list ladder: \n";
    */ /*
    std::cout << "vacuum Gamma:\n";
    glb_ainv = -2*sqrt(2.);
    glb_muc = 0.9;
    glb_mud = 1.;
    comp K1p_vac, K1a_vac, Gamma_T, Gamma_K1;
    double Gam0;
    K1p_vac = fRG_solve_K1r(0.,0.,'p',1e4,1e-8,3,0);
    K1a_vac = fRG_solve_K1r(0.,0.,'a',1e4,1e-8,3,0);
    Gamma_T = - 2*M_PI/(glb_mc*glb_md/(glb_mc+glb_md)*glb_ainv);
    Gam0 = -gint(1e4,1e-8,3);
    Gamma_K1 = Gam0 + K1p_vac + K1a_vac;
    //Gamma_K1 = fRG_solve_K1full(0.,0.,1e4,1e-8,1,0);
    std::cout << "K1p = " << K1p_vac << "\n";
    std::cout << "K1a = " << K1a_vac << "\n";
    std::cout << "Gamma_0 = " << Gam0 << "\n";
    //std::cout << "Gamma_K1 = " << Gamma_K1 << "\n";
    std::cout << "Gamma_T = " << Gamma_T << "\n";
    K1p_vac = ladder_K1r(0.,0.,'p',1e4,1e-8,1,0);
    K1a_vac = ladder_K1r(0.,0.,'a',1e4,1e-8,1,0);
    Gamma_T = - 2*M_PI/(glb_mc*glb_md/(glb_mc+glb_md)*glb_ainv);
    Gam0 = -gint(1e4,1e-8,1);
    Gamma_K1 = Gam0 + K1p_vac + K1a_vac;
    //Gamma_K1 = fRG_solve_K1full(0.,0.,1e4,1e-8,1,0);
    std::cout << "K1p = " << K1p_vac << "\n";
    std::cout << "K1a = " << K1a_vac << "\n";
    std::cout << "Gamma_0 = " << Gam0 << "\n";
    //std::cout << "Gamma_K1 = " << Gamma_K1 << "\n";
    std::cout << "Gamma_T = " << Gamma_T << "\n";
    glb_muc = 1.;
    */

    //ladder_list('f',1e4,1e-10,1,0,-2*sqrt(2),2*sqrt(2),-5,1e-10,100,1e-15,301);
    //ladder_list('p',1e4,1e-10,1,0,-2*sqrt(2),2*sqrt(2),-5,1e-10,100,1e-15,301);
    //ladder_list('a',1e4,1e-10,1,0,-2*sqrt(2),2*sqrt(2),-5,1e-10,100,1e-15,21);


    //ladder_p_list_wq(10.,10.,1e4,1e-10,6,6);


    //fRG_p_list_wq(10.,10.,1e4,1e-10,101,101);

        /*
        comp solvefRG1 = solve_1lfRG_nsc (0.2, 0.1, 1e1, 1e-4);
        std::cout << "fRG solution = " << solvefRG1 << "\n";

        comp test0001 = solve_1lfRG_nsc;
        std::cout << "fRG solution = " << test0001 << "\n";
        */

    //

    // 1-LOOP FRG WITH SOFT REGULATOR
    /*
    comp Gamm_fRG_exact_sharp = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-6, 1, 0);
    comp Gamm_fRG_num_sharp = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-6, 1, 1);
    comp Gamm_ladder_exact = ladder(0.1,0.1,'p',1e4,1e-6,1,0);
    comp Gamm_ladder_num = ladder(0.1,0.1,'p',1e4,1e-6,1,1);
    comp Gamm_fRG_exact_soft = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-6, 3, 0);
    //comp Gamm_fRG_num_soft = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-6, 3, 1);
    std::cout << "ladder: int = analytical, reg = sharp: " << Gamm_ladder_exact << "\n";
    std::cout << "ladder: int = numerical, reg = sharp: " << Gamm_ladder_num << "\n";
    std::cout << "fRG: int = analytical, reg = sharp: " << Gamm_fRG_exact_sharp << "\n";
    std::cout << "fRG: int = numerical, reg = sharp: " << Gamm_fRG_num_sharp << "\n";
    std::cout << "fRG: int = analytical, reg = soft: " << Gamm_fRG_exact_soft << "\n";
    //std::cout << "fRG: int = numerical, reg = soft: " << Gamm_fRG_num_soft << "\n";
    */
    /*
    comp integral_dlPi0_1 = perform_dL_Pi0_vpp_integral (1e4, 0.1, 0.1, 'c', 'd', 'p', 0);
    comp integral_dlPi0_2 = perform_dL_Pi0_vpp_integral (1, 0.1, 0.1, 'c', 'd', 'p', 0);
    comp integral_dlPi0_3 = perform_dL_Pi0_vpp_integral (1e-6, 0.1, 0.1, 'c', 'd', 'p', 0);
    comp integral_dlPi0_4 = perform_dL_Pi0_vpp_integral (1e-7, 0.1, 0.1, 'c', 'd', 'p', 0);

    std::cout << "soft v integral, Lambda = 1e4: " << integral_dlPi0_1 << "\n";
    std::cout << "soft v integral, Lambda = 1: " << integral_dlPi0_2 << "\n";
    std::cout << "soft v integral, Lambda = 1e-6: " << integral_dlPi0_3 << "\n";
    std::cout << "soft v integral, Lambda = 1e-7: " << integral_dlPi0_4 << "\n";
    //comp integral_dlPi0_5 = perform_dL_Pi0_vpp_integral (1e-10, 0.1, 0.1, 'c', 'd', 'p', 0);
    //std::cout << "soft v integral, Lambda = 0: " << integral_dlPi0_5 << "\n";
    */

    /*
    Test_soft_v_function<comp> test_soft_v_function(1.0,0.1,0,0);
    Test_soft_v_function<comp> test_soft_v_functionuoo(1.0,0.1,1,1);
    Test_soft_v_function<comp> test_soft_v_functionloo(1.0,0.1,-1,2);
    comp test_soft_v_f_integral = integrator<comp>(test_soft_v_function, -1, 1);
    comp test_soft_v_f_integraluoo = integrator<comp>(test_soft_v_functionuoo, 0, 1);
    comp test_soft_v_f_integralloo = integrator<comp>(test_soft_v_functionloo, 0, 1);
    comp test_soft_v_f_total = test_soft_v_f_integral+test_soft_v_f_integraluoo+test_soft_v_f_integralloo;
    std::cout << "Lambda = 1, w = 0.1: " << test_soft_v_f_total << "\n";

    comp function001 = functiontestsoft(-10.0, 1.0,0.1);
    std::cout << "function: " << function001 << "\n";
    */

    /*
    Integrand_dL_Pi0_vpp<comp> integrand_dL_Pi0_vpp13(10000,0.1,0.0,0.0,'c','d','p',0,0);
    comp integrand_test = integrand_dL_Pi0_vpp13(1.0);
    std::cout << "test integrand1: " << integrand_test << "\n";
    comp integral_test = perform_dL_Pi0_vpp_integral (10000, 0.1, 0.0, 'c', 'd', 'p', 0);
    std::cout << "test integral2: " << integral_test << "\n";
    integral_test = perform_dL_Pi0_vpp_integral (10000, 0.0, 0.0, 'c', 'd', 'p', 0);
    std::cout << "test integral3: " << integral_test << "\n";
    integral_test = perform_dL_Pi0_vpp_integral (1e-6, 0.0, 0.0, 'c', 'd', 'p', 0);
    std::cout << "test integral4: " << integral_test << "\n";
    */

    /*
    glb_mud = -8.0;
    glb_ainv = 0.0;
    Gamm_fRG_exact_soft = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-6, 3, 0);
    std::cout << "fRG: int = analytical, reg = soft: " << Gamm_fRG_exact_soft << "\n";

    Gamm_ladder_exact = ladder(0.1,0.1,'p',1e4,1e-6,1,0);
    std::cout << "ladder: int = analytical, reg = sharp: " << Gamm_ladder_exact << "\n";
    */

    double t1 = get_time();
    glb_ainv = sqrt(2.);
    glb_mud = -4.;
    //ladder_list('p',1e4,1e-6,1,0,-2.*sqrt(2),2.*sqrt(2),-10.,1e-10,100,1e-14,51);
    //ladder_list('p',1e4,1e-6,1,1,-2.*sqrt(2),2.*sqrt(2),-10.,1e-10,100,1e-14,51);
    //fRG_p_list(1e4,1e-6,1,0,-2.*sqrt(2),2.*sqrt(2),-10.,0.5,100,1e-14,51);
    //fRG_p_list(1e4,1e-6,1,1,-2.*sqrt(2),2.*sqrt(2),-10.,0.5,100,1e-14,51);
    //fRG_p_list(1e4,1e-6,3,0,-2.*sqrt(2),2.*sqrt(2),-10.,0.5,100,1e-14,51);
    //fRG_p_list(1e4,1e-6,3,1,-2.*sqrt(2),2.*sqrt(2),-10.,0.5,100,1e-14,27);
    //comp ladder_test = ladder(0.,0.1,'p',1e4,1e-6,1,1);
    //std::cout << "ladder: " << ladder_test << "\n";
    //comp Gamm_fRG_num_soft = fRG_solve_nsc(0., 0.1, 1e4, 1e-6, 3, 1);
    //std::cout << "fRG: int = numerical, reg = soft: " << Gamm_fRG_num_soft << "\n";
    //double find_root_fRG_test_reg3 = find_root_fRG (0., 0., 1e4, 1e-6, 3, 0, -3.0, sqrt(2.0), 0.5, 100, 1e-10);
    //std::cout << "mu_d = " << find_root_fRG_test_reg3 << "\n";
    get_time(t1);

    // A-CHANNEL
    // ======================================
    glb_ainv = -1.0;
    glb_mud = 0.0;
    comp gamma_p, gamma_a, Gamma;
    double Gamma0, dt;
    dt = get_time();

    // zero bosonic momentum
    /*
    std::cout << "ladder with analytical bubble integral:\n";
    gamma_p = ladder_K1r(0.,0.,'p',1e4,1e-8,1,0);
    gamma_a = ladder_K1r(0.,0.,'a',1e4,1e-8,1,0);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = ladder_full(0.,0.,1e4,1e-8,1,0);
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt); */
    /*
    comp bubble_intp = perform_Pi0_vpp_integral (0., 0., 'd', 'c', 'p', 1e4, 1e-8, 0);
    comp bubble_inta = perform_Pi0_vpp_integral (0., 0., 'd', 'c', 'a', 1e4, 1e-8, 0);
    std::cout << "bubble_a = " << bubble_inta << "\n";
    std::cout << "bubble_p = " << bubble_intp << "\n";
    comp exact_bubble_a = exact_bare_bubble (1., -3., 0., 'd', 'c', 'a');
    comp exact_bubble_p = exact_bare_bubble (1., -3., 0., 'd', 'c', 'p');
    std::cout << "exact bubble_a = " << exact_bubble_a << "\n";
    std::cout << "exact bubble_p = " << exact_bubble_p << "\n";
    */ /*
    dt = get_time();
    gamma_p = ladder_K1r(0.,0.,'p',1e4,1e-8,1,1);
    gamma_a = ladder_K1r(0.,0.,'a',1e4,1e-8,1,1);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = ladder_full(0.,0.,1e4,1e-8,1,1);
    std::cout << "ladder with numerical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.,'p',1e4,1e-8,1,0);
    gamma_a = fRG_solve_K1r(0.,0.,'a',1e4,1e-8,1,0);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = fRG_solve_K1full(0.,0.,1e4,1e-8,1,0);
    std::cout << "fRG with sharp regulator and analytical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.,'p',1e4,1e-8,1,1);
    gamma_a = fRG_solve_K1r(0.,0.,'a',1e4,1e-8,1,1);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = fRG_solve_K1full(0.,0.,1e4,1e-8,1,1);
    std::cout << "fRG with sharp regulator and numerical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.,'p',1e4,1e-8,3,0);
    gamma_a = fRG_solve_K1r(0.,0.,'a',1e4,1e-8,3,0);
    Gamma0 = -gint(1e4,1e-8,3);
    Gamma = fRG_solve_K1full(0.,0.,1e4,1e-8,3,0);
    std::cout << "fRG with soft regulator and analytical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.,'p',1e4,1e-8,3,1);
    gamma_a = fRG_solve_K1r(0.,0.,'a',1e4,1e-8,3,1);
    Gamma0 = -gint(1e4,1e-8,3);
    Gamma = fRG_solve_K1full(0.,0.,1e4,1e-8,3,1);
    std::cout << "fRG with soft regulator and analytical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt); */

    // finite bosonic momentum q != 0
    /*
    dt = get_time();
    std::cout << "ladder with analytical bubble integral:\n";
    gamma_p = ladder_K1r(0.,0.1,'p',1e4,1e-8,1,0);
    gamma_a = ladder_K1r(0.,0.1,'a',1e4,1e-8,1,0);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = ladder_full(0.,0.1,1e4,1e-8,1,0);
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = ladder_K1r(0.,0.1,'p',1e4,1e-8,1,1);
    gamma_a = ladder_K1r(0.,0.1,'a',1e4,1e-8,1,1);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = ladder_full(0.,0.1,1e4,1e-8,1,1);
    std::cout << "ladder with numerical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.1,'p',1e4,1e-8,1,0);
    gamma_a = fRG_solve_K1r(0.,0.1,'a',1e4,1e-8,1,0);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = fRG_solve_K1full(0.,0.1,1e4,1e-8,1,0);
    std::cout << "fRG with sharp regulator and analytical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.1,'p',1e4,1e-8,1,1);
    gamma_a = fRG_solve_K1r(0.,0.1,'a',1e4,1e-8,1,1);
    Gamma0 = -gint(1e4,1e-8,1);
    Gamma = fRG_solve_K1full(0.,0.1,1e4,1e-8,1,1);
    std::cout << "fRG with sharp regulator and numerical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);

    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.1,'p',1e4,1e-8,3,0);
    gamma_a = fRG_solve_K1r(0.,0.1,'a',1e4,1e-8,3,0);
    Gamma0 = -gint(1e4,1e-8,3);
    Gamma = fRG_solve_K1full(0.,0.1,1e4,1e-8,3,0);
    std::cout << "fRG with soft regulator and analytical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);
    */ /*
    dt = get_time();
    gamma_p = fRG_solve_K1r(0.,0.1,'p',1e4,1e-8,3,1);
    gamma_a = fRG_solve_K1r(0.,0.1,'a',1e4,1e-8,3,1);
    Gamma0 = -gint(1e4,1e-8,3);
    Gamma = fRG_solve_K1full(0.,0.1,1e4,1e-8,3,1);
    std::cout << "fRG with soft regulator and numerical bubble integral:\n";
    std::cout << "K1_a = " << gamma_a << "\n";
    std::cout << "K1_p = " << gamma_p << "\n";
    std::cout << "Gamma0 = " << Gamma0 << "\n";
    std::cout << "K1 = " << Gamma << "\n";
    get_time(dt);
    */
    /*
    double gint_test;
    glb_sharpness = 2.0;
    gint_test = gint(1e4, 1e-10, 1);
    std::cout << "reg = 1, g = " << gint_test << "\n";
    gint_test = gint(1e4, 1e-10, 2);
    std::cout << "reg = 2, g = " << gint_test << "\n";
    gint_test = gint(1e4, 1e-10, 3);
    std::cout << "reg = 3, g = " << gint_test << "\n";
    gint_test = gint(1e4, 1e-10, 4);
    std::cout << "reg = 4, g = " << gint_test << "\n";
    gint_test = gint(1e4, 1e-10, 5);
    std::cout << "reg = 5, g = " << gint_test << "\n";
    gint_test = gint(1e4, 1e-10, 6);
    std::cout << "reg = 6, g = " << gint_test << "\n";

    double test33 = std::tgamma(0.25)*(4.-pow(2.,3./4.))/(8*M_PI*M_PI);
    std::cout << "test = " << test33 << "\n";
    */
    // TEST GORKOV (CONST. GAMMA)

    glb_muc = 1.0;
    glb_mud = -10.0;
    glb_ainv = -1.0;
    // comp gorkov_test = fRG_solve_K1r(0., sqrt(2.*glb_muc), 'c', 1e4, 1e-10, 1, 0);
    // std::cout << "Gamma = " << gorkov_test << "\n";

    //std::vector<double> myvec{3.12,3.23,42.564,43.23,567.75,23.54};
    //std::vector<double>::iterator test_begin = myvec.begin();
    //std::cout<< "begin = " << myvec.begin() << "\n";

    // TEST PAID-INTEGRATOR
    // =================================s
    /*
    std::cout << "\nTest the PAID-integrator: \n";
    Domain1D<comp> am(0.,0.5);
    Domain1D<comp> mb(0.5, 1.);
    PAIDInput test_integrand_paid1 (am, f_testpaid, 0);
    PAIDInput test_integrand_paid2 (mb, f_testpaid, 0);

    //f_testpard(3.0+2.0);
    std::vector<PAIDInput> test_integrands_paid = {test_integrand_paid1, test_integrand_paid2};

    PAID test_integral_paid(test_integrands_paid);


    std::cout << "integral test: \n";
    auto paid_solution = test_integral_paid.solve()[0];
    std::cout << "integral test: \n";
    std::cout << "integral = " << paid_solution << "\n";

    Domain1D<comp> ab01(0.,1.);
    PAIDInput integrand_paid_gauss (ab01,integrand_semiinfinite_gauss,0);
    PAID integral_paid_gauss({integrand_paid_gauss});
    std::cout << "gauss-integral = " << integral_paid_gauss.solve()[0] << "\n";

    std::pair<int, double> pair_test = {1, 2.3};
    std::cout << pair_test.first  << " and " << pair_test.second << "\n";
    std::map<int, double> m { {0, 1.0}, {1,exp(1.)},{2,exp(2.)}, {3,exp(3.)}, };
    std::cout << "m(3) = " << m[3] << "\n";

    Domain1D<comp> abm11(-1.,1.);
    Integrand_Pi0_theta<comp> integrandPi0Theta_paid(0.5,-0.1,0.3,0.4,'c','d');
    PAIDInput integrandPi0Theta_paid_input(abm11,integrandPi0Theta_paid,0);
    PAID integralPi0Theta_paid({integrandPi0Theta_paid_input});
    std::cout << "theta-integral paid = " << 1./(4.*M_PI*M_PI)*0.4*0.4*integralPi0Theta_paid.solve()[0] << "\n";
    comp keldysh_theta_integral_result = perform_integral_Pi0_theta (0.5,-0.1,0.3,0.4,'c','d',1);
    std::cout << "theta-integral gauss-lobatto = " << keldysh_theta_integral_result << "\n";
    */

    paid::Domain<1> d({0.},{1.});
    double a = 0.1;
    double exact_result = 1./sqrt(a);
    Integrand_Gauss gauss(a);
    paid::PAIDInput<1,Integrand_Gauss, int> paid_integrand{d,gauss,0};
    paid::PAIDConfig config;
    paid::PAID<1,Integrand_Gauss,double,int,double> paid_integral(config);
    std::complex<double> paid_result = paid_integral.solve({paid_integrand})[0];
    std::cout << "Gauss integral for a = " << a << ": exact: 1/a^0.5 = " << exact_result << ", paid =" << paid_result << "\n";

    paid::Domain<2> d2({-M_PI, -M_PI},{M_PI, M_PI});
    Integrand_sin2D sinsin(1.);
    paid::PAIDInput<2,Integrand_sin2D,int> paid_integrand2d{d2,sinsin,0};
    paid::PAIDConfig config2;
    paid::PAID<2,Integrand_sin2D,double,int,std::array<double,2>> paid_integral2d(config2);
    std::complex<double> paid_result2d = paid_integral2d.solve({paid_integrand2d})[0];
    std::cout << "2D sin integral: " << paid_result2d << "\n";


    /*
    std::vector<double> v_test(10,0.0);
    std::cout << "v_test = (" << v_test[0];
    for (int i = 1; i<v_test.size(); i++){
        std::cout << ", " << v_test[i];
    }
    std::cout << ")\n";
    glb_mud = 0.0;
    double wmax = 10., vppmax = 10., qmax = 10., kmax = 10., vpp, kpp, t_kpp;
    int nw = 6, nvpp = 6, nq = 3, nk = 51;
    comp exact, lobatto, paid;
    double dt_inttype;
    dt_inttype = get_time();
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2*vppi*vppmax/(nw-1);
            for (int qi = 0; qi < nq; ++qi) {
                q = qi*qmax/(nq-1);
                /*
                for (int kppi = 0; kppi < nk; ++kppi){
                    kpp = kppi*kmax/(nk-1);
                    lobatto = perform_integral_Pi0_theta(w,vpp,q,kpp,'c','d',1);
                    std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", k = " << kpp << ", lobatto = " << lobatto << "\n";
                    paid = perform_integral_Pi0_theta(w,vpp,q,kpp,'c','d',2);
                    std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", k = " << kpp << ", paid = " << paid << "\n";
                }*/ /*

                for (int ti = 0; ti < nk; ++ti){
                    t_kpp = ti/(nk-1.);
                    Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo_lobatto(w, vpp, q, 0.0, 'c', 'd',1,1);
                    //Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo_paid(w, vpp, q, 0.0, 'c', 'd',2,1);

                    lobatto = integrand_Pi0_kpp_0oo_lobatto(t_kpp);
                    std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", t = " << t_kpp << ", lobatto = " << lobatto << "\n";
                    //paid = integrand_Pi0_kpp_0oo_paid(t_kpp);
                    //std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", t = " << t_kpp << ", paid = " << paid << "\n";
                }

                exact = exact_bare_bubble(w, vpp, q, 'c', 'd','a');
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", exact = " << exact << "\n";
                //lobatto = perform_integral_Pi0_kpp_chan(w, vpp, q, 'c', 'd', 1, 'a');
                //std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", lobatto = " << lobatto << "\n";
                //paid = perform_integral_Pi0_kpp_chan(w,vpp,q,'c','d',2,'a');
                //std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", paid = " << paid << "\n";
            }
        }
    } //*/ /*
    get_time(dt_inttype);
    dt_inttype = get_time();
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2*vppi*vppmax/(nw-1);
            for (int qi = 0; qi < nq; ++qi) {
                q = qi*qmax/(nq-1);

                for (int ti = 0; ti < nk; ++ti){
                    t_kpp = ti/(nk-1.);
                    //Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo_lobatto(w, vpp, q, 0.0, 'c', 'd',1,1);
                    Integrand_Pi0_kpp<comp> integrand_Pi0_kpp_0oo_paid(w, vpp, q, 0.0, 'c', 'd',2,1);

                    //lobatto = integrand_Pi0_kpp_0oo_lobatto(t_kpp);
                    //std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", t = " << t_kpp << ", lobatto = " << lobatto << "\n";
                    paid = integrand_Pi0_kpp_0oo_paid(t_kpp);
                    std::cout << "v1 = " << w << ", v2 = " << vpp << ", q = " << q << ", t = " << t_kpp << ", paid = " << paid << "\n";
                }

                exact = exact_bare_bubble(w, vpp, q, 'c', 'd','a');
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", exact = " << exact << "\n";
                //lobatto = perform_integral_Pi0_kpp_chan(w, vpp, q, 'c', 'd', 1, 'a');
                //std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", lobatto = " << lobatto << "\n";
                //paid = perform_integral_Pi0_kpp_chan(w,vpp,q,'c','d',2,'a');
                //std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", paid = " << paid << "\n";
            }
        }
    } //*/ /*
    get_time(dt_inttype);
    dt_inttype = get_time();
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2 * wi * wmax / (nw - 1);
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2 * vppi * vppmax / (nw - 1);
            for (int qi = 0; qi < nq; ++qi) {
                q = qi * qmax / (nq - 1);
                lobatto = perform_integral_Pi0_kpp_chan(w, vpp, q, 'c', 'd', 1, 'a');
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", lobatto = " << lobatto << "\n";
            }
        }
    }
    get_time(dt_inttype);
    dt_inttype = get_time();
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2 * wi * wmax / (nw - 1);
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2 * vppi * vppmax / (nw - 1);
            for (int qi = 0; qi < nq; ++qi) {
                q = qi * qmax / (nq - 1);
                paid = perform_integral_Pi0_kpp_chan(w,vpp,q,'c','d',2,'a');
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", paid = " << paid << "\n";
            }
        }
    }
    get_time(dt_inttype);
    */
    std::vector<double> v(8,1.0);
    std::cout << "v1 = (" << v[0];
    for (int i=1; i<v.size(); ++i){
        std::cout << ", " << v[i];
    }
    std::cout << ")\n";

    std::vector<double> v2{5,21,3,1,32,43,2,51,6,20,35,6};
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";
    std::cout << "after make_heap:\n";
    std::make_heap(v2.begin(),v2.end());
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";
    std::pop_heap(v2.begin(),v2.end());
    std::cout << "after pop_heap:\n";
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";
    v2.pop_back();
    std::cout << "after pop_back:\n";
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";
    v2.push_back(34);
    std::cout << "push_back:\n";
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";
    std::push_heap(v2.begin(), v2.end());
    std::cout << "push_heap:\n";
    std::cout << "v2 = (" << v2[0];
    for (int i=1; i<v2.size(); ++i){
        std::cout << ", " << v2[i] ;
    }
    std::cout << ")\n";


    //list_bubble_int_tkpp(10.,1.,1,-5.,1.,'c','d',101,11,11,6);
    //list_bubble_int_tkpp(10.,1.,2,-5.,1.,'c','d',101,11,11,6);

    // INTEGRATE BARE BUBBLE NUMERICALLY
    // =================================

    //list_bubble_int_theta (-2.0, 1.0, 0.01, 0.1, 'c', 'd', 100);
    //list_bubble_int_kpp (10.,10.,-5,1,'c','d',100,13,13,13);

    /*
    comp theta_integral_001 = perform_integral_Pi0_theta (0.0, 0.0, 0., 0.1, 'c', 'd');
    std::cout << "theta-integral = " << theta_integral_001 << "\n";
    comp theta_integral_002 = perform_integral_Pi0_theta (0.0, 0.0, 1e-10, 0.1, 'c', 'd');
    std::cout << "theta-integral = " << theta_integral_002 << "\n";
    comp theta_integral_003 = perform_integral_Pi0_theta (0.0, 0.0, -1e-10, 0.1, 'c', 'd');
    std::cout << "theta-integral = " << theta_integral_003 << "\n";
    comp bubble_integral_001 = perform_integral_Pi0_kppt(0.1,0.1,0.1,'c','d');
    comp bubble_exact_001 = exact_bare_bubble_v1v2(0.1,0.1,0.1,'c','d');
    std::cout << "numerical: " << bubble_integral_001 << ", exact: " << bubble_exact_001 << "\n";
    comp bubble_integral_002 = perform_integral_Pi0_kppt(1e-4,0.0,0.0,'c','d');
    comp bubble_exact_002 = exact_bare_bubble_v1v2(1e-4,0.0,0.0,'c','d');
    std::cout << "numerical: " << bubble_integral_002 << ", exact: " << bubble_exact_002 << "\n";
    comp bubble_integral_003 = perform_integral_Pi0_kppt(-1e-4,0.,0.,'c','d');
    comp bubble_exact_003 = exact_bare_bubble_v1v2(-1e-4,0.,0.,'c','d');
    std::cout << "numerical: " << bubble_integral_003 << ", exact: " << bubble_exact_003 << "\n";
    comp bubble_exact_004 = exact_bare_bubble_v1v2(0.,0.,0.,'c','d');
    comp bubble_integral_004 = perform_integral_Pi0_kppt(0.0,0.0,0.,'c','d');
    std::cout << "numerical: " << bubble_integral_004 << ", exact: " << bubble_exact_004 << "\n";

    std::cout << "now test \n";
    comp test_integral_Pi0_001 = perform_integral_Pi0_kpp_chan (-9.6, 4.8, 10, 'c', 'd','p');
    std::cout <<"test result = " << test_integral_Pi0_001 << "\n";
    test_integral_Pi0_001 = perform_integral_Pi0_kpp_chan (-9.6, 4.8, 0, 'c', 'd','p');
    std::cout <<"test result = " << test_integral_Pi0_001 << "\n";
    test_integral_Pi0_001 = perform_integral_Pi0_kpp_chan (-0.2, 0.1, 9, 'c', 'd','p');
    std::cout <<"test result = " << test_integral_Pi0_001 << "\n";


    double wmax = 10.;
    double vppmax = 10.;
    double qmax = 10.;
    int nw = 201;
    int nvpp = 201;
    int nq = 11;
    double vpp;
    comp result_integral;
    */
    /*
    for (int wi = 0; wi < nw; ++wi) {
        w = -wmax + 2*wi*wmax/(nw-1);
        for (int vppi = 0; vppi < nvpp; ++vppi) {
            vpp = -vppmax + 2*vppi*vppmax/(nw-1);
            for (int qi = 0; qi < nq; ++qi) {
                q = qi*qmax/(nq-1);
                if ((std::abs(w/2 + vpp) < 1e-15) || (std::abs(w/2 - vpp) < 1e-15)) {
                    std::cout << "v1 and/or v2 = 0 \n";
                    result_integral = 0.0;
                }
                else {
                    result_integral = perform_integral_Pi0_kpp_chan (w, vpp, q, 'c', 'd', 'p');

                }
                std::cout << "w = " << w << ", vpp = " << vpp << ", q = " << q << ", result = " << result_integral << "\n";
            }
        }
    }
     */

    std::cout << "exact bare bubble : \n";
    //integral_bubble_w_vpp_list_exact('c','d','p',10.,10.,10.,201,201,6);
    std::cout << "numerical bare bubble : \n";
    //integral_bubble_w_vpp_list_integrator('c','d','p',10.0,10.0,10.0,201,201,6);

    // test Gauss as integral with infinite intervals
    /*comp result_gauss = integrator<comp>(integrand_infinite_gauss, 0.0, 1.0);
    std::cout << "Gauss integral = " << result_gauss << "\n";
    result_gauss = integrator<comp>(integrand_semiinfinite_gauss, 0.0, 1.0);
    std::cout << "Gauss integral = " << result_gauss << "\n"; */

    // integrate bare bubble with regulator
    //comp test_regulated_integral = perform_integral_Pi0_vpp (1.0, 0.0, 0.0, 'c', 'd', 'p');
    //std::cout << "regulated bubble = " << test_regulated_integral << "\n";

    // TEST INTEGRATOR AND SOLVER
    /*
    std::cout << "test integrate \n";


    comp integrate_result;
    integrate_result = perform_SimpleBubble_integral (0.03, -2.0, 0.2, 0.4, 'c', 'c');
    std::cout << "x-integral = " << integrate_result << "\n";

    comp y_fin;
    const comp y_ini = SimpleBubble(0.03, -2.0, 0.2, 0.4, -1.0, 'c', 'c');
    const double x_ini = -1.0;
    const double x_fin = 1.0;

    ODE_solver_RK4(y_fin, x_ini, y_ini, x_fin,
                   rhs_test3, sq_substitution, sq_resubstitution, nODE);



    std::cout << "yfin = " << y_fin << "\n";
    std::cout << "nODE = " << nODE << "\n";
    */

    // LOOP INTEGRAL SELFENERGY
    glb_mud = -2.5;
    glb_ainv = sqrt(2.);
    /*
    comp testloop001 = selfenergy_ladder (0.1, 0.0,1e4,1e-10);
    std::cout << "loop selfenergy test = " << testloop001 << "\n";
    comp testloop002 = selfenergy_ladder (0.0, 0.0,1e4,1e-10);
    std::cout << "loop selfenergy test = " << testloop002 << "\n";
    */
    //selfenergy_ladder_list_vk(10.,10.,1e4,1e-10,11,6);

    // Try out things concerning PAID integrator

    std::array<std::size_t,4> i_list{126,0,127, 64};
    std::size_t composite_index_i_list;
    composite_index_i_list = paid::get_composite_index<4>(128,i_list);
    std::cout << "composite_index = " << composite_index_i_list << "\n";
    std::array<std::size_t,4> i_list_recovered;
    i_list_recovered = paid::get_each_index<4>(128,composite_index_i_list);
    std::cout << "i_list_recovered = (";
    for (std::size_t i = 0; i < size(i_list); ++i){
        std::cout << i_list_recovered[i] << ",";
    }
    std::cout << ")\n";

    std::array<std::size_t,1> i_list1d{126};
    std::size_t composite_index_i_list1d;
    composite_index_i_list1d = paid::get_composite_index<1>(128,i_list1d);
    std::cout << "composite_index = " << composite_index_i_list1d << "\n";
    std::array<std::size_t,1> i_list_recovered1d;
    i_list_recovered1d = paid::get_each_index<1>(128,composite_index_i_list1d);
    std::cout << "i_list_recovered = (";
    for (std::size_t i = 0; i < size(i_list1d); ++i){
        std::cout << i_list_recovered1d[i] << ",";
    }
    std::cout << "\n";

    std::array<std::size_t,3> k_vector{0,0,0};
    for (std::size_t k = 0; k < 8; ++k) {
        k_vector = paid::get_each_index<3>(2,k);
        std::cout << "k = << " << k << ", k_vector = (" << k_vector[0] << ", " << k_vector[1] << ", " << k_vector[2] << ")\n";
    }

    std::cout << ")\n";

    get_time(t0);

    std::cout << "Goodbye World! \n";

    if (MPI_FLAG) {
        MPI_Finalize();
    }

}



