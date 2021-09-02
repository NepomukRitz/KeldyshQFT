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
#include "Monte-Carlo_Trial.h"
#include "fRG-T-matrix-approach.h"
#include "1D-integrals.h"
#include "../utilities/util.h"
#include "../ODE_solvers.h"
#include "../grids/flow_grid.h"



int main() {

    if (MPI_FLAG) {
        MPI_Init(nullptr, nullptr);
    }

    gsl_set_error_handler_off();

    std::cout << "Hello world!\n";

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
    glb_prec = 1e-16;

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
    std::cout << "Bubble integral = " << bubble1 << ", exact = " << bubble1e << /* ", simplified = " << output2t << */  "\n";
    //double test1 = 1e4;
    //double test2 = sqrt(test1);
    //std::cout << test1 << " = " << test2 << "\n";
    std::cout << "Bubble integral = " << bubble2 << ", exact = " << bubble2e << /* ", simplified = " << output3t << */ "\n";
    std::cout << "Bubble integral = " << bubble3 << ", exact = " << bubble3e << /* ", simplified = " << output4t << */ "\n";
    std::cout << "Bubble integral = " << bubble4 << ", exact = " << bubble4e << /* ", simplified = " << output5t << */ "\n";
    std::cout << "test now \n";
    glb_muc = 0.0;
    glb_mud = 0.0;
    glb_ainv = 1.0;

    comp ladder1, ladder2, ladder3, ladder4, ladder5;
    double laddere;
    std::cout << "k integral exact\n";
    ladder1 = ladder(0.0,0.0,10,1e-10,1,0);
    ladder2 = ladder(0.0,0.0,1e4,1e-10,1,0);
    ladder3 = ladder(0.0,0.0,10,1e-5,1,0);
    ladder4 = ladder(0.0,0.0,1e2,1e-5,1,0);
    ladder5 = ladder(0.0,0.0,1e4,1e-5,1,0);
    laddere = -4*M_PI;
    std::cout << "ladder = " << ladder1 << "\n";
    std::cout << "ladder = " << ladder2 << "\n";
    std::cout << "ladder = " << ladder3 << "\n";
    std::cout << "ladder = " << ladder4 << "\n";
    std::cout << "ladder = " << ladder5 << "\n";
    std::cout << "exact = " << laddere << "\n";



    std::cout << "k integral numerical\n";
    ladder1 = ladder(0.0,0.0,10,1e-10,1,1);
    ladder2 = ladder(0.0,0.0,1e4,1e-10,1,1);
    ladder3 = ladder(0.0,0.0,10,1e-5,1,1);
    ladder4 = ladder(0.0,0.0,1e2,1e-5,1,1);
    ladder5 = ladder(0.0,0.0,1e4,1e-5,1,1);
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
    comp ladder01 = 1./ladder(0.0,0.0,1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    glb_mud = -1.23;
    ladder01 = 1./ladder(0.0,0.0,1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    glb_mud = -1.2325;
    ladder01 = 1./ladder(0.0,0.0,1e4,1e-10,1,0);
    std::cout << "a^(-1) = " << glb_ainv << ", mu_d = " << glb_mud << ": " << ladder01 << "\n";
    double root_test0 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,-2.,1e-10,100,1e-16);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test0 << "\n";


    glb_ainv = 2.0;
    double root_test1 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,-2.,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test1 << "\n";
    glb_ainv = 1.5;
    double root_test2 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test1,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test2 << "\n";
    glb_ainv = 1.0;
    double root_test3 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test2,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test3 << "\n";
    glb_ainv = 0.5;
    double root_test4 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test3,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test4 << "\n";
    glb_ainv = 0.0;
    double root_test5 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test4,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test5 << "\n";
    glb_ainv = -0.5;
    double root_test6 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test5,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test6 << "\n";
    glb_ainv = -1.0;
    double root_test7 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test6,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test7 << "\n";
    glb_ainv = -2.0;
    double root_test8 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,root_test7,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test8 << "\n";
    glb_ainv = 0.282843;
    double root_test9 = find_root_newton(0.0,0.0,1e4,1e-10,1,0,-0.342722,1e-10,100,1e-10);
    std::cout << "a^(-1) = " << glb_ainv << ": " << root_test9 << "\n";

    std::cout << "now solver with exact momentum integral:\n";
    ladder_p_list (1e4, 1e-10, 1, 0, -2.*sqrt(2), 2.*sqrt(2), -10., 1e-10, 100, 1e-16, 301);
    std::cout << "now solver with numerical momentum integral:\n";
    //ladder_p_list (1e4,1e-10,1,1,-2*sqrt(2),2*sqrt(2),-10.,1e-10,100,1e-16,101);


    // comp testfmu = test_f_mu(0.0,0.0,1e3,1e-4,1,1.0);
    // std::cout << "f_mu = " << testfmu << "\n";

    // 1-LOOP FRG
    // ==================================

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

    comp Gamma_fRG_exact = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-10, 1, 0);
    comp Gamma_fRG_num = fRG_solve_nsc(0.1, 0.1, 1e4, 1e-10, 1, 1);
    comp Gamma_ladder_exact = ladder(0.1,0.1,1e4,1e-10,1,0);
    comp Gamma_ladder_num = ladder(0.1,0.1,1e4,1e-10,1,1);

    std::cout << "Gamma_fRG exact k-integral: " << Gamma_fRG_exact << "\n";
    std::cout << "Gamma_fRG numerical k-integral: " << Gamma_fRG_num << "\n";
    std::cout << "Gamma_ladder exact k-integral: " << Gamma_ladder_exact << "\n";
    std::cout << "Gamma_ladder numerical k-integral: " << Gamma_ladder_num << "\n";

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
    fRG_p_list (1e4, 1e-10, 1,1, -2*sqrt(2), 2*sqrt(2), -10.0, 0.5, 100, 1e-16, 11);

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


    std::cout << "now list fRG: \n";
    //fRG_p_list_wq(10.,10.,1e4,1e-10,201,101);
    std::cout << "now list ladder: \n";
    //ladder_p_list_wq(10.,10.,1e4,1e-10,6,6);


    //fRG_p_list_wq(10.,10.,1e4,1e-10,101,101);

        /*
        comp solvefRG1 = solve_1lfRG_nsc (0.2, 0.1, 1e1, 1e-4);
        std::cout << "fRG solution = " << solvefRG1 << "\n";

        comp test0001 = solve_1lfRG_nsc;
        std::cout << "fRG solution = " << test0001 << "\n";
        */
    //

    // INTEGRATE BARE BUBBLE NUMERICALLY
    // =================================


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
    double result_gauss = integrator<double>(integrand_infinite_gauss, 0.0, 1.0);
    std::cout << "Gauss integral = " << result_gauss << "\n";
    result_gauss = integrator<double>(integrand_semiinfinite_gauss, 0.0, 1.0);
    std::cout << "Gauss integral = " << result_gauss << "\n";

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


    get_time(t0);

    std::cout << "Goodbye World! \n";

    if (MPI_FLAG) {
        MPI_Finalize();
    }

}



