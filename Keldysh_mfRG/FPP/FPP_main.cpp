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
#include "1D-integrals.h"
#include "../utilities/util.h"
#include "../ODE_solvers.h"



int main() {

#ifdef MPI_FLAG
    MPI_Init(nullptr, nullptr);
#endif

    gsl_set_error_handler_off();

    cout << "Hello world!\n";

    double t0 = get_time();

    // Coordinate Transformations
    // =======================================

    /*double x;
    double y;
    cout << "Type in x.\n";
    cin >> x;
    cout << "Type in y.\n";
    cin >> y;
    double angle;
    angle = phi(x, y);
    double angledegree = angle*180/M_PI;
    cout << "Phi is " << angledegree << ".\n";
    return 0;*/
    //coordinate_transform();
    //coordinate_transform();
    /*vector<double> v_in;
    v_in = {1.0,1.0,1.0};
    vector<double> v_out;
    v_out = cartesion_to_spherical(v_in);
    cout << v_out << ".\n"; */
    /* for (int n = 1; n < 100; ++n) {
        vector<double> v_in = {1.1223498, -2.24342338, 3.322433};
        vector<double> v_out;
        v_out = transform_back_forth(v_in, n);
        double vx = v_out[0];
        double vy = v_out[1];
        double vz = v_out[2];
        cout << "N = " << n << ": v = {" << setprecision(8) << vx << "," << setprecision(8) << vy << "," << setprecision(8) << vz << "}.\n";
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

    // Cout Bubble
    // =======================================

    /*cout << "Give nu, k^2, i\n";
    double nu;
    double ksquared;
    char i;
    cin >> nu;
    cin >> ksquared;
    cin >> i;
    cout << "Your bare Green's function is :" << G0(nu, ksquared, i) << "\n";
    cout << "Give omega, nu, k, q, cos(theta), channel, i, j\n";
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
    cout << "Your bare bubble is: " << SimpleBubble(omega, nu, q, k, x, channel, i, j) << "\n"; */

    // Calculate Bubble integrals
    // =======================================

    //integrate_bubble_full_monte_carlo(0.0,0.0,0.0,'c','d','p',1000.,50000,0.1);
    //comp output = integrate_bubble_vegas(0.0,0.0,0.0,'c','d','p',1000.,50000,0.1);
    //cout << "Result: " << output << "\n";

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

        cout << "result: " << res << ", error: " << err << "\n";
    }
    //*/

    // Test if clause with char
    //====================================
    /*
    char c;
    cout << "Type in character: ";
    cin >> c;
    double output2 = test2(c);
    cout << "\n Your value is " << output2 << ".\n";
    */

    // 1D INTEGRALS
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
    cout << test001 << "\n";
    comp test002 = pow(-3.0+glb_i,0.5);
    cout << test002 << "\n";
    comp test003 = -1./(2*sqrt(2)*M_PI*pow(-1.0+1e-16*glb_i,0.5));
    cout << test003 << "\n";
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
    cout << "Output 1: " << output1 << "\n";
    output2 = bubble_integrate_kpp (0.1, 0.1, 0.1, 'c', 'c', 0, 0.,0.);
    cout << "Output 2: " << output2 << "\n";
    output3 = bubble_k_2d_integrated (0.1, 0.1, 0.1, 'c', 'c', 'p', 0, 0);
    double reoutput3, imoutput3;
    reoutput3 = real(output3);
    imoutput3 = imag(output3);
    cout << "Output 3: " << reoutput3 << " + i " << imoutput3 << "\n";
    print_numerical_bubble(0.1, 0.1, 0.1, 'c', 'c', 'p', 0, 0);*/


    //cout << "now new" << "\n";

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
    cout << test_heaviside << "\n";
    test_heaviside = heaviside(x_neg);
    cout << test_heaviside << "\n";
    test_heaviside = heaviside(x_zero);
    cout << test_heaviside << "\n";*/


    /*comp test_sfbb_a, test_sfbb_p, test_sfbb_t;
    test_sfbb_a = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 'a');
    cout << test_sfbb_a << "\n";
    test_sfbb_p = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 'p');
    cout << test_sfbb_p << "\n";
    test_sfbb_t = sharp_frequency_exact_bare_bubble(-2.3, 0.4, 1.2,'c', 'd', 't');
    cout << test_sfbb_t << "\n";

    cout << "now solving ODE: \n";
    double w = 3.4, q = 1.2, g = -0.5, Lambda_i = 100.0, Lambda_f = 0.001;
    comp outputK1p = K1cdcd_solution(w, q, g, Lambda_i, 'p',Lambda_f, 1e-10, 1e-10, 0.0);
    comp outputK1a = K1cdcd_solution(w, q, g, Lambda_i, 'a',Lambda_f, 1e-10, 1e-10, 0.0);
    comp outputK1t = K1cdcd_solution(w, q, g, Lambda_i, 't',Lambda_f, 1e-10, 1e-10, 0.0);
    cout << "K1p = " << outputK1p << "\n";
    cout << "K1a = " << outputK1a << "\n";
    cout << "K1t = " << outputK1t << "\n";
    w = 0.0, q = 0.0, g = -0.1, Lambda_i = 100.0, Lambda_f = 1e-10;
    outputK1p = K1cdcd_solution(w, q, g, Lambda_i, 'p',Lambda_f, 1e-10, 1e-10, 0.0);
    outputK1a = K1cdcd_solution(w, q, g, Lambda_i, 'a',Lambda_f, 1e-10, 1e-10, 0.0);
    outputK1t = K1cdcd_solution(w, q, g, Lambda_i, 't',Lambda_f, 1e-10, 1e-10, 0.0);
    cout << "K1p = " << outputK1p << "\n";
    cout << "K1a = " << outputK1a << "\n";
    cout << "K1t = " << outputK1t << "\n";
    *//*outputK1p = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1a = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1t = K1cdcd_solution_nint(0.1, 0.1, -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    cout << "K1p = " << outputK1p << "\n";
    cout << "K1a = " << outputK1a << "\n";
    cout << "K1t = " << outputK1t << "\n";*/
    /*
    cout << "now exact with w = 0 = q: \n";
    outputK1p = K1cdcd_solution(0., 0., -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0);
    outputK1a = K1cdcd_solution(0., 0., -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0);
    outputK1t = K1cdcd_solution(0., 0., -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0);
    cout << "K1p = " << outputK1p << "\n";
    cout << "K1a = " << outputK1a << "\n";
    cout << "K1t = " << outputK1t << "\n";
    //K1Lambda (0.0, 0.0, -0.1, 'p', 10.0, 1e-10, 100, 1e-12, 1e-12, 0.0);
    //K1Lambdag(0.0,0.0,-0.1,0.1,'p',10.0,1e-10,100,1e-12,1e-12,0.0);
    *//*outputK1p = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 'p',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1a = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 'a',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    outputK1t = K1cdcd_solution_nint(0., 0., -0.1, 10.0, 't',0.1, 1e-10, 1e-10, 0.0, 1e-12,1e-12);
    cout << "K1p = " << outputK1p << "\n";
    cout << "K1a = " << outputK1a << "\n";
    cout << "K1t = " << outputK1t << "\n";*/

    //K1Lambda (0.1, 0.1, -0.1, 'p', 10.0, 1e-10, 100, 1e-12, 1e-12, 0.0, 1e-12, 1e-12);

        /*derivativefsqrt();
        ODE_solver_test();
        double outputre, outputim;
        comp outputcomp;
        outputre = dLsfebb ( 0.1, 1.0, 0.1, 'c', 'c', 'p', 0, 1e-10);
        outputim = dLsfebb ( 0.1, 1.0, 0.1, 'c', 'c', 'p', 1, 1e-10);
        outputcomp = outputre + glb_i * outputim;
        cout << "Derivative is " << outputcomp << "\n";
        ODE_solver_K1p(3.0, 3.0, 0.1, 10.0, -0.1, 10);*/

    // Sharp frequency regulator
    // =============================

    /*comp output13 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'a');
    comp output14 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'p');
    comp output15 = sharp_frequency_exact_bare_bubble( 2., 3., 4., 'c','d' ,'t');
    cout << "Sharp regulator: " << output13 << ".\n";
    cout << "Sharp regulator: " << output14 << ".\n";
    cout << "Sharp regulator: " << output15 << ".\n";*/


    // DOUBLE PENDULUM
    // ===============================

    /*cout << "now double pendulum: " << "\n";
    solve_double_pendulum(9.81,1.0,2.0,1.0,1.5,0.0,10.0,0.0,30.0,1.0,-2.3);*/

    // ZEROS
    /*
    double testroot1 = find_root_divergence (sigmoidal, -1.0, 1., 1000, 10e-10);
    cout << "x = " << testroot1 << "\n";
    comp exactvalue = sigmoidal(testroot1);
    cout << "sigmoid(x) = " << exactvalue << "\n";

    double testroot2 = fine_root_newton (fzeros1, 3, 10e-10, 1000, 10e-10);
    cout << "x = " << testroot2 << "\n";
    double exactvalue2 = fzeros1(testroot2);
    cout << "f1(x) = " << exactvalue2 << "\n";

    double testroot3 = fine_root_newton (fzeros2, -1.0, 1., 1000, 10e-10);
    cout << "x = " << testroot3 << "\n";
    double xx0 = pow(M_PI,3)/log(2)/(10*M_EULER);
    cout << "x_exact " << xx0 << "\n";
    double exactvalue3 = fzeros2(testroot3);
    cout << "f2(x) = " << exactvalue3 << "\n";
    */

    // LADDER APPROXIMATION

    glb_muc = 0.0;
    glb_mud = 0.0;
    glb_ainv = 1.0;

    comp bubble1 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 10, 1e-10);
    double bubble1e = exactzerobubble(10, 1e-10);
    //comp output2t = perform_vacuum_integral (10, 1e-10);
    comp bubble2 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 1e4, 1e-10);
    double bubble2e = exactzerobubble(1e4, 1e-10);
    //comp output3t = perform_vacuum_integral (1, 1, 1e4, 1e-10);
    comp bubble3 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 10, 1e-5);
    double bubble3e = exactzerobubble(10, 1e-5);
    //comp output4t = perform_vacuum_integral (1, 1, 10, 1e-5);
    comp bubble4 = perform_Pi0_vpp_integral (0.0, 0.0, 'c', 'd', 'p', 1e4, 1e-5);
    double bubble4e = exactzerobubble(1e4, 1e-5);
    //comp output5t = perform_vacuum_integral (1, 1, 1e4, 1e-5); */
    cout << "Bubble integral = " << bubble1 << ", exact = " << bubble1e << /* ", simplified = " << output2t << */  "\n";
    //double test1 = 1e4;
    //double test2 = sqrt(test1);
    //cout << test1 << " = " << test2 << "\n";
    cout << "Bubble integral = " << bubble2 << ", exact = " << bubble2e << /* ", simplified = " << output3t << */ "\n";
    cout << "Bubble integral = " << bubble3 << ", exact = " << bubble3e << /* ", simplified = " << output4t << */ "\n";
    cout << "Bubble integral = " << bubble4 << ", exact = " << bubble4e << /* ", simplified = " << output5t << */ "\n";

    glb_muc = 0.0;
    glb_mud = 0.0;
    glb_ainv = 1.0;

    comp ladder1, ladder2, ladder3, ladder4, ladder5;
    double laddere;
    ladder1 = ladder(0.0,0.0,10,1e-10,1);
    ladder2 = ladder(0.0,0.0,1e4,1e-10,1);
    ladder3 = ladder(0.0,0.0,10,1e-5,1);
    ladder4 = ladder(0.0,0.0,1e2,1e-5,1);
    ladder5 = ladder(0.0,0.0,1e4,1e-5,1);
    laddere = -4*M_PI;
    cout << "ladder = " << ladder1 << "\n";
    cout << "ladder = " << ladder2 << "\n";
    cout << "ladder = " << ladder3 << "\n";
    cout << "ladder = " << ladder4 << "\n";
    cout << "ladder = " << ladder5 << "\n";
    cout << "exact = " << laddere << "\n";

    // comp testfmu = test_f_mu(0.0,0.0,1e3,1e-4,1,1.0);
    // cout << "f_mu = " << testfmu << "\n";

    // TEST INTEGRATOR AND SOLVER
    /*
    cout << "test integrate \n";


    comp integrate_result;
    integrate_result = perform_SimpleBubble_integral (0.03, -2.0, 0.2, 0.4, 'c', 'c');
    cout << "x-integral = " << integrate_result << "\n";



    comp y_fin;
    const comp y_ini = SimpleBubble(0.03, -2.0, 0.2, 0.4, -1.0, 'c', 'c');
    const double x_ini = -1.0;
    const double x_fin = 1.0;

    ODE_solver_RK4(y_fin, x_ini, y_ini, x_fin,
                   rhs_test, sq_substitution, sq_resubstitution, nODE);



    cout << "yfin = " << y_fin << "\n";
    cout << "nODE = " << nODE << "\n";
     */

    get_time(t0);

    cout << "Goodbye World! \n";

#ifdef MPI_FLAG
    MPI_Finalize();
#endif

}



