#ifndef KELDYSH_MFRG_TESTFUNCTIONS_H
#define KELDYSH_MFRG_TESTFUNCTIONS_H

#include <cmath>                    // use M_PI as pi
#include "../propagator.h"                // propagators
#include "../state.h"                  // State class
#include "../loop.h"                   // self-energy loop
#include "../bubbles.h"                // bubble function
#include "../ODE_solvers.h"               // ODE solvers
#include "../right_hand_sides.h"       // compute the right hand sides of flow equations
#include "../utilities/write_data2file.h"        // writing data to txt or hdf5 file
#include "../utilities/hdf5_routines.h"          // writing states to hdf5 file
#include "../perturbation_theory.h"
#include <boost/math/special_functions/polygamma.hpp> // Polygamma function
#include "../integrator/integrator.h"
#include "../postprocessing/causality_FDT_checks.h"   // check causality and FDTs




/**
 * Function to test the loop function and the calculation of the SelfEnergy
 * @param state : State initialized with initial conditions
 */
template <typename Q>
void testSelfEnergy_and_K1(State<Q>& state, double Lambda){

    Propagator<Q> g(Lambda, state.selfenergy, 'g');

    //Calculate the vertex
    sopt_state(state, Lambda);

    Vertex<Q> temp_vertex_a (1, Lambda), temp_vertex_p (1, Lambda); //All zeros
    temp_vertex_a[0].avertex() = state.vertex[0].avertex();
    temp_vertex_p[0].pvertex() = state.vertex[0].pvertex();

    loop(state.selfenergy, temp_vertex_a, g, false);//Calculate the SelfEnergy in SOPT

    //Print results in .h5 format
    cvec SER(nFER);
    cvec SEK(nFER);

    for(int iv = 0; iv<nFER; iv++){
        SER[iv] = state.selfenergy.val(0, iv, 0);
        SEK[iv] = state.selfenergy.val(1, iv, 0);
    }

    //Print results
    write_h5_rvecs("SOPT_test.h5", {"v", "ReSER", "ImSER", "ReSEK", "ImSEK",
                                  "K1a_R", "K1a_I", "K1p_R", "K1p_I", "K1t_R", "K1t_I"},
                                 {ffreqs, SER.real(), SER.imag(), SEK.real(), SEK.imag(),
                                  state.vertex[0].avertex().K1.real(),
                                  state.vertex[0].avertex().K1.imag(),
                                  state.vertex[0].pvertex().K1.real(),
                                  state.vertex[0].pvertex().K1.imag(),
                                  state.vertex[0].tvertex().K1.real(),
                                  state.vertex[0].tvertex().K1.imag()});

}


/**
 * Function that prints out a .h5 file with the value of the rhs of a SOPT flow at the given Lambda for both an FFT and an fRG calculation
 * @param Lambda    : Lambda at which the derivatives are to be calculated
 */
template <typename Q>
void test_derivatives_K1a(double Lambda){
    vec<Q> blah(nw1_a);
    vec<Q> rhs_SOPT_FFT_K1a = dSOPT_FFT_K1a_rhs(blah, Lambda);
    vec<Q> rhs_flow = rhs_bubbles_flow(blah, Lambda);

    write_h5_rvecs("derivatives_K1a.h5",
                   {"v", "FFT_R", "FFT_I", "SOPT_R", "SOPT_I"},
                   {ffreqs, rhs_SOPT_FFT_K1a.real(), rhs_SOPT_FFT_K1a.imag(), rhs_flow.real(), rhs_flow.imag()});
}

#if MAX_DIAG_CLASS >= 2
/**
 * This function checks the consistency of the K2 class
 * @param Lambda    : Scale
 * @param r         : Channel to be tested i.e. K_2^r
 */
template <typename Q>
auto test_K2_consistency(double Lambda, const char r) -> bool{
    State<Q> bare (Lambda);   //Create a bare state
    bare.initialize();  //Initialize bare state

    Propagator<Q> G(Lambda, bare.selfenergy, 'g'); //Bare propagator at scale Lambda>

    //Create and initialize states to save the channel contributions in
    State<Q> K1a (Lambda);
    State<Q> K1p (Lambda);
    State<Q> K1t (Lambda);
    K1a.initialize();
    K1p.initialize();
    K1t.initialize();

    //Calculate K1a, K1p and K1t contributions separately
    bubble_function(K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false);

    //Create and initialize the K2r-objects to test
    State<Q> test_K2r_with_K1a (Lambda);
    State<Q> test_K2r_with_K1p (Lambda);
    State<Q> test_K2r_with_K1t (Lambda);
    test_K2r_with_K1a.initialize();
    test_K2r_with_K1p.initialize();
    test_K2r_with_K1t.initialize();


    if(r=='p' ||  r=='t') {
        //Perform TOPT calculation of K2r with K1a
        bubble_function(test_K2r_with_K1a.vertex, K1a.vertex, bare.vertex, G, G, r, false);
    }

    if(r=='a' || r=='t') {
        //Perform TOPT calculation of K2r with K1p
        bubble_function(test_K2r_with_K1p.vertex, K1p.vertex, bare.vertex, G, G, r, false);
    }

    if(r=='a' || r=='p') {
        //Perform TOPT calculation of K2r with K1t
        bubble_function(test_K2r_with_K1t.vertex, K1t.vertex, bare.vertex, G, G, r, false);
    }


    //TOPT-consistency checks for K2

    bool empty = true;  //Boolean that stays true if everything is zero
    if(r=='a' || r=='p') {
#pragma omp parallel
        //Check parallelized that everything in the K2 vertex is zero.
        for (int index = 0; index < test_K2r_with_K1t.vertex[0].avertex().K2.size(); ++index) {
            if (test_K2r_with_K1t.vertex[0].avertex().K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) print("TOPT-consistency check passed. K2a or K2p with K1t is zero everywhere.", true);
        else print("TOPT-consistency check failed. K2a or K2p with K1t is not zero everywhere.", true);
    }
    else if(r=='t'){
#pragma omp parallel
        //Check parallelized that everything in the K2 vertices is zero.
        for (int index = 0; index < test_K2r_with_K1a.vertex[0].avertex().K2.size(); ++index) {
            if (test_K2r_with_K1p.vertex[0].tvertex().K2_acc(index) != 0.) {
                empty = false;  //If any one element is not zero, K2 is not empty
                break;  //Exit the for-loop early
            }
        }
        //Print result of consistency check
        if (empty) print("TOPT-consistency check passed. K2t with K1a and K1p are zero everywhere.", true);
        else print("TOPT-consistency check failed. K2t with K1a and K1p are not zero everywhere.", true);
    }

    return empty;
}
#endif

// to check central part of multi-loop flow equations:
// compute diagrams with non-symmetric intermediate results
void compute_non_symmetric_diags(const double Lambda, bool write_flag = false) {
    State<state_datatype> bare (Lambda); // bare state
    bare.initialize();         // initialize bare state

    Propagator<state_datatype> G (Lambda, bare.selfenergy, 'g'); // bare propagator
    Propagator<state_datatype> S (Lambda, bare.selfenergy, 's'); // bare differentiated propagator = single scale propagator

    // Psi := K1p in PT2 + bare vertex
    State<state_datatype> Psi (Lambda);
    Psi.initialize();         // initialize bare state
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, G, G, 'p', false);

    // K1a_dot in PT2
    State<state_datatype> PT2_K1adot (Lambda);
    bubble_function(PT2_K1adot.vertex, bare.vertex, bare.vertex, G, S, 'a', true);
    // K1p_dot in PT2
    State<state_datatype> PT2_K1pdot (Lambda);
    bubble_function(PT2_K1pdot.vertex, bare.vertex, bare.vertex, G, S, 'p', true);

    if (write_flag) {
        write_hdf("Psi_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, Psi);
        write_hdf("PT2_K1a_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1adot);
        write_hdf("PT2_K1p_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1pdot);
    }

    std::vector<State<state_datatype>> central_bubblestates = {PT2_K1adot, PT2_K1pdot};

    for (int i = 0; i < 2; i++){
        State<state_datatype> centralstate_dot = central_bubblestates[i];

        // intermediate results
        State<state_datatype> K1rdot_PIa_K1p (Lambda);
        bubble_function(K1rdot_PIa_K1p.vertex, centralstate_dot.vertex, Psi.vertex, G, G, 'a', false);

        State<state_datatype> K1p_PIa_K1rdot (Lambda);
        bubble_function(K1p_PIa_K1rdot.vertex, Psi.vertex, centralstate_dot.vertex, G, G, 'a', false);


        if (write_flag) {
            write_hdf("K1rdot_PIa_K1p_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_hdf("K1p_PIa_K1rdot_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
        }

        Vertex<state_datatype> dGammaL_half1 = K1rdot_PIa_K1p.vertex;
        Vertex<state_datatype> dGammaR_half1 = K1p_PIa_K1rdot.vertex;
        dGammaL_half1[0].half1().reorder_due2antisymmetry(dGammaR_half1[0].half1());
        dGammaR_half1[0].half1().reorder_due2antisymmetry(dGammaL_half1[0].half1());
        K1rdot_PIa_K1p.vertex = dGammaL_half1;
        K1p_PIa_K1rdot.vertex = dGammaR_half1;

        // create non-symmetric vertex with differentiated vertex on the left
        GeneralVertex<state_datatype , non_symmetric> dGammaL(n_spin, Lambda);
        dGammaL[0].half1()  = dGammaL_half1[0].half1();  // assign half 1 to dGammaL
        dGammaL[0].half2() = dGammaR_half1[0].half1();  // assign half 2 as half 1 to dGammaR [symmetric -> left()=right()]

        // insert this non-symmetric vertex on the right of the bubble
        State<state_datatype> dGammaC_r(Lambda);
        bubble_function(dGammaC_r.vertex, Psi.vertex, dGammaL, G, G, 'a', false);


        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<state_datatype , non_symmetric> dGammaR (n_spin, Lambda);
        dGammaR[0].half1() = dGammaR_half1[0].half1();  // assign half 1
        dGammaR[0].half2() = dGammaL_half1[0].half1();  // assign half 2 as half 1 of dGammaL

        // insert this non-symmetric vertex on the left of the bubble
        State<state_datatype> dGammaC_l(Lambda);
        bubble_function(dGammaC_l.vertex, dGammaR, Psi.vertex, G, G, 'a', false);


        if (write_flag) {
            write_hdf("K1rdot_PIa_K1p_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_hdf("K1p_PIa_K1rdot_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
            write_hdf("dGammaC_r_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_r);
            write_hdf("dGammaC_l_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_l);
        }
    }
}

/**
 * Function that computes K1 (and K2, K3) up to PT4, and performs FDT checks
 */
void test_PT4(double Lambda, bool write_flag = false) {
    print("Compute K1 (and K2, K3) up to PT4.", true);
    // Initialize a bare state
    State<state_datatype> bare (Lambda);
    bare.initialize();

    // Initialize a bare propagator
    Propagator<state_datatype> G(Lambda, bare.selfenergy, 'g');

    // Compute K1 in PT2
    State<state_datatype> PT2_K1a (Lambda);
    State<state_datatype> PT2_K1p (Lambda);
    State<state_datatype> PT2_K1t (Lambda);

    double t0 = get_time();
    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT2_K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false);
    print("Computed K1 in PT2.", true);
    get_time(t0);
    if (write_flag) {
        write_hdf("PT2_K1a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1a);
        write_hdf("PT2_K1p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1p);
        write_hdf("PT2_K1t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1t);
    }

    // Compute K1 in PT3, using K1 in PT2
    State<state_datatype> PT3_K1a (Lambda);
    State<state_datatype> PT3_K1p (Lambda);
    State<state_datatype> PT3_K1t (Lambda);

    t0 = get_time();
    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT3_K1p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'p', false);
    // for K1t in PT3, need a-vertex in PT2 due to a <-> t symmetry
    bubble_function(PT3_K1t.vertex, PT2_K1t.vertex + PT2_K1a.vertex, bare.vertex, G, G, 't', false);
#if MAX_DIAG_CLASS >= 2
    // set K2 part of this vertex to zero
    PT3_K1t.vertex[0].tvertex().K2 = vec<state_datatype> (PT3_K1t.vertex[0].tvertex().K2.size());
#endif
    print("Computed K1 in PT3.", true);
    get_time(t0);
    if (write_flag) {
        write_hdf("PT3_K1a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1a);
        write_hdf("PT3_K1p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1p);
        write_hdf("PT3_K1t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K1t);
    }

    // Compute K2 in PT3, using K1p, K1t in PT2
    State<state_datatype> PT3_K2a (Lambda);
    State<state_datatype> PT3_K2p (Lambda);
    State<state_datatype> PT3_K2t (Lambda);
    State<state_datatype> PT3_K2t_a (Lambda);
    State<state_datatype> PT3_K2t_p (Lambda);

    t0 = get_time();
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'a', false);   // K2a in PT3 (PT2_K1t = 0)
    bubble_function(PT3_K2p.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'p', false);   // K2p in PT3 (PT2_K1t = 0)
    bubble_function(PT3_K2t_a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 't', false);   // contribution of K2t in PT3 obtained by inserting K1a in PT2
    bubble_function(PT3_K2t_p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 't', false);   // contribution of K2t in PT3 obtained by inserting K1p in PT2
    // PT3_K2t_a should also have a K1-contribution due to a-t symmetry (PT2_K1t implicitly inserted) --> set to zero
    PT3_K2t_a.vertex[0].tvertex().K1 = vec<state_datatype> (PT3_K2t_a.vertex[0].tvertex().K1.size());
    PT3_K2t.vertex = PT3_K2t_a.vertex + PT3_K2t_p.vertex; // sum of contributions from a- and p-insertions

    // K2' in PT3 would be obtained by flipping the left and right vertex, but since K2' is not saved, these terms would give zero
    print("Computed K2 in PT3.", true);
    get_time(t0);
    if (write_flag) {
        write_hdf("PT3_K2a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2a);
        write_hdf("PT3_K2p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2p);
        write_hdf("PT3_K2t_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t);
        write_hdf("PT3_K2t_a_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t_a);
        write_hdf("PT3_K2t_p_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT3_K2t_p);
    }


    // Compute K3 in PT4, using K1 and K2 in PT2 and PT3
    State<state_datatype> PT4_31 (Lambda);
    State<state_datatype> PT4_31_a_a1 (Lambda);
    State<state_datatype> PT4_31_a_p1 (Lambda);
    State<state_datatype> PT4_31_a_t1 (Lambda);
    State<state_datatype> PT4_31_a_a2 (Lambda);
    State<state_datatype> PT4_31_a_p2 (Lambda);
    State<state_datatype> PT4_31_a_t2 (Lambda);
    State<state_datatype> PT4_31_p_a1 (Lambda);
    State<state_datatype> PT4_31_p_p1 (Lambda);
    State<state_datatype> PT4_31_p_t1 (Lambda);
    State<state_datatype> PT4_31_p_a2 (Lambda);
    State<state_datatype> PT4_31_p_p2 (Lambda);
    State<state_datatype> PT4_31_p_t2 (Lambda);
    State<state_datatype> PT4_31_t_a1 (Lambda);
    State<state_datatype> PT4_31_t_p1 (Lambda);
    State<state_datatype> PT4_31_t_t1 (Lambda);
    State<state_datatype> PT4_31_t_a2 (Lambda);
    State<state_datatype> PT4_31_t_p2 (Lambda);
    State<state_datatype> PT4_31_t_t2 (Lambda);

    State<state_datatype> PT4_13 (Lambda);
    State<state_datatype> PT4_13_a_a1 (Lambda);
    State<state_datatype> PT4_13_a_p1 (Lambda);
    State<state_datatype> PT4_13_a_t1 (Lambda);
    State<state_datatype> PT4_13_a_a2 (Lambda);
    State<state_datatype> PT4_13_a_p2 (Lambda);
    State<state_datatype> PT4_13_a_t2 (Lambda);
    State<state_datatype> PT4_13_p_a1 (Lambda);
    State<state_datatype> PT4_13_p_p1 (Lambda);
    State<state_datatype> PT4_13_p_t1 (Lambda);
    State<state_datatype> PT4_13_p_a2 (Lambda);
    State<state_datatype> PT4_13_p_p2 (Lambda);
    State<state_datatype> PT4_13_p_t2 (Lambda);
    State<state_datatype> PT4_13_t_a1 (Lambda);
    State<state_datatype> PT4_13_t_p1 (Lambda);
    State<state_datatype> PT4_13_t_t1 (Lambda);
    State<state_datatype> PT4_13_t_a2 (Lambda);
    State<state_datatype> PT4_13_t_p2 (Lambda);
    State<state_datatype> PT4_13_t_t2 (Lambda);

    State<state_datatype> PT4_22 (Lambda);
    State<state_datatype> PT4_22_a_aa (Lambda);
    State<state_datatype> PT4_22_a_ap (Lambda);
    State<state_datatype> PT4_22_a_pa (Lambda);
    State<state_datatype> PT4_22_a_pp (Lambda);
    State<state_datatype> PT4_22_p_aa (Lambda);
    State<state_datatype> PT4_22_p_ap (Lambda);
    State<state_datatype> PT4_22_p_pa (Lambda);
    State<state_datatype> PT4_22_p_pp (Lambda);
    State<state_datatype> PT4_22_t_aa (Lambda);
    State<state_datatype> PT4_22_t_ap (Lambda);
    State<state_datatype> PT4_22_t_pa (Lambda);
    State<state_datatype> PT4_22_t_pp (Lambda);

    t0 = get_time();

    // a-channel:
    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_a_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT4_31_a_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT4_31_a_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT4_31_a_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT4_31_a_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT4_31_a_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 'a', false);

    if (write_flag) {
        write_hdf("PT4_31_a_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_a1);
        write_hdf("PT4_31_a_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_p1);
        write_hdf("PT4_31_a_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_t1);
        write_hdf("PT4_31_a_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_a2);
        write_hdf("PT4_31_a_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_p2);
        write_hdf("PT4_31_a_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_a_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex[0].avertex() = PT4_31_a_a1.vertex[0].avertex()
                             + PT4_31_a_p1.vertex[0].avertex()
                             + PT4_31_a_t1.vertex[0].avertex()
                             + PT4_31_a_a2.vertex[0].avertex()
                             + PT4_31_a_p2.vertex[0].avertex()
                             + PT4_31_a_t2.vertex[0].avertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    // this can only give a K1 contribution, since we do not compute K2'
    bubble_function(PT4_13_a_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false);
    // the following should all give zero
    bubble_function(PT4_13_a_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'a', false);
    bubble_function(PT4_13_a_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 'a', false);
    bubble_function(PT4_13_a_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false);
    bubble_function(PT4_13_a_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 'a', false);
    bubble_function(PT4_13_a_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 'a', false);

    if (write_flag) {
        write_hdf("PT4_13_a_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_a1);
        write_hdf("PT4_13_a_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_p1);
        write_hdf("PT4_13_a_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_t1);
        write_hdf("PT4_13_a_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_a2);
        write_hdf("PT4_13_a_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_p2);
        write_hdf("PT4_13_a_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_a_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex[0].avertex() = PT4_13_a_a1.vertex[0].avertex()
                             + PT4_13_a_p1.vertex[0].avertex()
                             + PT4_13_a_t1.vertex[0].avertex()
                             + PT4_13_a_a2.vertex[0].avertex()
                             + PT4_13_a_p2.vertex[0].avertex()
                             + PT4_13_a_t2.vertex[0].avertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_a_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false);
    bubble_function(PT4_22_a_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 'a', false);
    bubble_function(PT4_22_a_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 'a', false);
    bubble_function(PT4_22_a_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'a', false);

    if (write_flag) {
        write_hdf("PT4_22_a_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_aa);
        write_hdf("PT4_22_a_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_ap);
        write_hdf("PT4_22_a_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_pa);
        write_hdf("PT4_22_a_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_a_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex[0].avertex() = PT4_22_a_aa.vertex[0].avertex()
                             + PT4_22_a_ap.vertex[0].avertex()
                             + PT4_22_a_pa.vertex[0].avertex()
                             + PT4_22_a_pp.vertex[0].avertex();

    print("Computed a-channel in PT4.", true);

    // p-channel:
    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_p_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT4_31_p_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT4_31_p_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT4_31_p_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT4_31_p_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT4_31_p_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 'p', false);

    if (write_flag) {
        write_hdf("PT4_31_p_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_a1);
        write_hdf("PT4_31_p_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_p1);
        write_hdf("PT4_31_p_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_t1);
        write_hdf("PT4_31_p_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_a2);
        write_hdf("PT4_31_p_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_p2);
        write_hdf("PT4_31_p_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_p_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex[0].pvertex() = PT4_31_p_a1.vertex[0].pvertex()
                             + PT4_31_p_p1.vertex[0].pvertex()
                             + PT4_31_p_t1.vertex[0].pvertex()
                             + PT4_31_p_a2.vertex[0].pvertex()
                             + PT4_31_p_p2.vertex[0].pvertex()
                             + PT4_31_p_t2.vertex[0].pvertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    // this can only give a K1 contribution, since we do not compute K2'
    bubble_function(PT4_13_p_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'p', false);
    // the following should all give zero
    bubble_function(PT4_13_p_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'p', false);
    bubble_function(PT4_13_p_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 'p', false);
    bubble_function(PT4_13_p_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'p', false);
    bubble_function(PT4_13_p_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 'p', false);
    bubble_function(PT4_13_p_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 'p', false);

    if (write_flag) {
        write_hdf("PT4_13_p_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_a1);
        write_hdf("PT4_13_p_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_p1);
        write_hdf("PT4_13_p_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_t1);
        write_hdf("PT4_13_p_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_a2);
        write_hdf("PT4_13_p_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_p2);
        write_hdf("PT4_13_p_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_p_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex[0].pvertex() = PT4_13_p_a1.vertex[0].pvertex()
                             + PT4_13_p_p1.vertex[0].pvertex()
                             + PT4_13_p_t1.vertex[0].pvertex()
                             + PT4_13_p_a2.vertex[0].pvertex()
                             + PT4_13_p_p2.vertex[0].pvertex()
                             + PT4_13_p_t2.vertex[0].pvertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_p_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'p', false);
    bubble_function(PT4_22_p_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 'p', false);
    bubble_function(PT4_22_p_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 'p', false);
    bubble_function(PT4_22_p_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'p', false);

    if (write_flag) {
        write_hdf("PT4_22_p_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_aa);
        write_hdf("PT4_22_p_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_ap);
        write_hdf("PT4_22_p_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_pa);
        write_hdf("PT4_22_p_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_p_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex[0].pvertex() = PT4_22_p_aa.vertex[0].pvertex()
                             + PT4_22_p_ap.vertex[0].pvertex()
                             + PT4_22_p_pa.vertex[0].pvertex()
                             + PT4_22_p_pp.vertex[0].pvertex();

    print("Computed p-channel in PT4.", true);

    // t-channel:
    // in the t-channel, we need to insert a and t simultaneously due to a <-> t symmetry // TODO: remove?
    // (spin sum in the t-channel makes use of this symmetry)                             // TODO: remove?

    // 3-1: insert all possible PT3 vertices on the left, bare vertex on the right
    bubble_function(PT4_31_t_a1.vertex, PT3_K1a.vertex, bare.vertex, G, G, 't', false);
    bubble_function(PT4_31_t_p1.vertex, PT3_K1p.vertex, bare.vertex, G, G, 't', false);
    bubble_function(PT4_31_t_t1.vertex, PT3_K1t.vertex, bare.vertex, G, G, 't', false);
    bubble_function(PT4_31_t_a2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 't', false);
    bubble_function(PT4_31_t_p2.vertex, PT3_K2p.vertex, bare.vertex, G, G, 't', false);
    bubble_function(PT4_31_t_t2.vertex, PT3_K2t.vertex, bare.vertex, G, G, 't', false);

    if (write_flag) {
        write_hdf("PT4_31_t_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_a1);
        write_hdf("PT4_31_t_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_p1);
        write_hdf("PT4_31_t_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_t1);
        write_hdf("PT4_31_t_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_a2);
        write_hdf("PT4_31_t_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_p2);
        write_hdf("PT4_31_t_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31_t_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_31.vertex[0].tvertex() = PT4_31_t_a1.vertex[0].tvertex()
                             + PT4_31_t_p1.vertex[0].tvertex()
                             + PT4_31_t_t1.vertex[0].tvertex()
                             + PT4_31_t_a2.vertex[0].tvertex()
                             + PT4_31_t_p2.vertex[0].tvertex()
                             + PT4_31_t_t2.vertex[0].tvertex();

    // 1-3: insert bare vertex on the left, all possible PT3 vertices on the right
    bubble_function(PT4_13_t_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 't', false);
    bubble_function(PT4_13_t_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 't', false);
    bubble_function(PT4_13_t_t1.vertex, bare.vertex, PT3_K1t.vertex, G, G, 't', false);
    bubble_function(PT4_13_t_a2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 't', false);
    bubble_function(PT4_13_t_p2.vertex, bare.vertex, PT3_K2p.vertex, G, G, 't', false);
    bubble_function(PT4_13_t_t2.vertex, bare.vertex, PT3_K2t.vertex, G, G, 't', false);

    if (write_flag) {
        write_hdf("PT4_13_t_a1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_a1);
        write_hdf("PT4_13_t_p1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_p1);
        write_hdf("PT4_13_t_t1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_t1);
        write_hdf("PT4_13_t_a2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_a2);
        write_hdf("PT4_13_t_p2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_p2);
        write_hdf("PT4_13_t_t2_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13_t_t2);
    }

    // sum of contributions obtained from different insertions
    PT4_13.vertex[0].tvertex() = PT4_13_t_a1.vertex[0].tvertex()
                             + PT4_13_t_p1.vertex[0].tvertex()
                             + PT4_13_t_t1.vertex[0].tvertex()
                             + PT4_13_t_a2.vertex[0].tvertex()
                             + PT4_13_t_p2.vertex[0].tvertex()
                             + PT4_13_t_t2.vertex[0].tvertex();

    // 2-2: insert all possible PT2 vertices on the left and right (PT2_K1t is always zero -> neglect it)
    bubble_function(PT4_22_t_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 't', false);
    bubble_function(PT4_22_t_ap.vertex, PT2_K1a.vertex, PT2_K1p.vertex, G, G, 't', false);
    bubble_function(PT4_22_t_pa.vertex, PT2_K1p.vertex, PT2_K1a.vertex, G, G, 't', false);
    bubble_function(PT4_22_t_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 't', false);

    if (write_flag) {
        write_hdf("PT4_22_t_aa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_aa);
        write_hdf("PT4_22_t_ap_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_ap);
        write_hdf("PT4_22_t_pa_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_pa);
        write_hdf("PT4_22_t_pp_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22_t_pp);
    }

    // sum of contributions obtained from different insertions
    PT4_22.vertex[0].tvertex() = PT4_22_t_aa.vertex[0].tvertex()
                             + PT4_22_t_ap.vertex[0].tvertex()
                             + PT4_22_t_pa.vertex[0].tvertex()
                             + PT4_22_t_pp.vertex[0].tvertex();


    print("Computed t-channel in PT4.", true);

    if (write_flag) {
        write_hdf("PT4_31_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_31);
        write_hdf("PT4_13_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_13);
        write_hdf("PT4_22_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT4_22);
    }

    print("Computed K1, K2, K3 in PT4.", true);
    get_time(t0);

    /** Make automated checks of all diagrams: Compute the values of all diagrams at all frequencies equal to zero,
     * and for all pairs of diagrams that should cancel, compute the relative deviation of their sum from zero.
     * Print all results to log.
     * */

    print("--- CHECK RESULTS: ---", true);
    print("--- print relative error of quantities that should be zero ---", true);

#ifdef KELDYSH_FORMALISM
    int iK = 1;
#else
    int iK = 0;
#endif
    // input variables: all frequencies equal to zero
    VertexInput input_a (iK, 0., 0., 0., 0, 0, 'a');
    VertexInput input_p (iK, 0., 0., 0., 0, 0, 'p');
    VertexInput input_t (iK, 0., 0., 0., 0, 0, 't');

#ifdef KELDYSH_FORMALISM
    vec<int> iK2s = {1, 2, 4}; // Keldysh indices of fully retarded components of K2
#else
    vec<int> iK2s = {0}; // Keldysh indices of Matsubara component of K2
#endif

    // labels to be printed to log
    std::string check_labels[] {"PT2: K1a + K1p: ", "PT2: K1t: ",
                           "PT2: K1a - exact: ", "PT2: K1p - exact: ",
                           "PT3: K1a - exact: ", "PT3: K1p - exact: ", "PT3: K1t - exact: ",
                           "PT3: K2a[1] - exact: ", "PT3: K2p[1] - exact: ", "PT3: K2t[1] - exact: ",
#ifdef KELDYSH_FORMALISM
                           "PT3: K2a[2] - exact: ", "PT3: K2p[2] - exact: ", "PT3: K2t[2] - exact: ",
                           "PT3: K2a[4] - exact: ", "PT3: K2p[4] - exact: ", "PT3: K2t[4] - exact: "
                            ,
#endif
                           "PT4: (K2a[1] <- K1p) + (K2p[1] <- K1a): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K1p) + (K2p[2] <- K1a): ",
                           "PT4: (K2a[4] <- K1p) + (K2p[4] <- K1a): ",
#endif
                           "PT4: (K2a[1] <- K1t) + (K2p[1] <- K1t): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K1t) + (K2p[2] <- K1t): ",
                           "PT4: (K2a[4] <- K1t) + (K2p[4] <- K1t): ",
#endif
                           "PT4: (K2a[1] <- K2a) + (K2p[1] <- K2p): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2a) + (K2p[2] <- K2p): ",
                           "PT4: (K2a[4] <- K2a) + (K2p[4] <- K2p): ",
#endif
                           "PT4: (K2a[1] <- K2p) + (K2p[1] <- K2a): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2p) + (K2p[2] <- K2a): ",
                           "PT4: (K2a[4] <- K2p) + (K2p[4] <- K2a): ",
#endif
                           "PT4: (K2a[1] <- K2t) + (K2p[1] <- K2t): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2a[2] <- K2t) + (K2p[2] <- K2t): ",
                           "PT4: (K2a[4] <- K2t) + (K2p[4] <- K2t): ",
#endif
                           "PT4: (K2t[1] <- K2a) + (K2t[1] <- K2t): ",
#ifdef KELDYSH_FORMALISM
                           "PT4: (K2t[2] <- K2a) + (K2t[2] <- K2t): ",
                           "PT4: (K2t[4] <- K2a) + (K2t[4] <- K2t): ",
#endif
                           "PT4: K3a + K3p: ",
                           "PT4: K3t (aa) + K3t (ap) + K3t (pa): "
                           };

    // K1 in PT2
    state_datatype PT2_K1a_0 = PT2_K1a.vertex[0].avertex().valsmooth<k1>(input_a, PT2_K1a.vertex[0].tvertex());
    state_datatype PT2_K1p_0 = PT2_K1p.vertex[0].pvertex().valsmooth<k1>(input_p, PT2_K1p.vertex[0].pvertex());
    state_datatype PT2_K1t_0 = PT2_K1t.vertex[0].tvertex().valsmooth<k1>(input_t, PT2_K1t.vertex[0].avertex());

    // K1 in PT3
    state_datatype PT3_K1a_0 = PT3_K1a.vertex[0].avertex().valsmooth<k1>(input_a, PT3_K1a.vertex[0].tvertex());
    state_datatype PT3_K1p_0 = PT3_K1p.vertex[0].pvertex().valsmooth<k1>(input_p, PT3_K1p.vertex[0].pvertex());
    state_datatype PT3_K1t_0 = PT3_K1t.vertex[0].tvertex().valsmooth<k1>(input_t, PT3_K1t.vertex[0].avertex());
#ifdef KELDYSH_FORMALISM
    state_datatype PT2_K1_exact = -(1./2.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 1);
    state_datatype PT3_K1_exact = -(1./2.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
#ifdef ZERO_TEMP
    state_datatype PT2_K1_exact = - glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 1);
    state_datatype PT3_K1_exact = - glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
    state_datatype PT2_K1_exact = - glb_U * pow(glb_U * boost::math::polygamma(1, (glb_Gamma + Lambda) / 2. / (glb_T * 2.*M_PI) + 0.5) / (glb_T * 2 * M_PI * M_PI), 1);
    state_datatype PT3_K1_exact = - glb_U * pow(glb_U * boost::math::polygamma(1, (glb_Gamma + Lambda) / 2. / (glb_T * 2.*M_PI) + 0.5) / (glb_T * 2 * M_PI * M_PI), 2);
#endif
#endif
    std::cout << "PT2 K1 exact: " << PT2_K1_exact << "\n";
    std::cout << "Computed value: " << PT2_K1a_0 << "\n";

    // K1 in PT4
    state_datatype PT4_K1a_0_ladder = PT4_31_a_a1.vertex[0].avertex().valsmooth<k1>(input_a, PT4_31_a_a1.vertex[0].tvertex());
    state_datatype PT4_K1p_0_ladder = PT4_31_p_p1.vertex[0].pvertex().valsmooth<k1>(input_p, PT4_31_p_p1.vertex[0].pvertex());
    state_datatype PT4_K1a_0_nonladder = PT4_13_a_a2.vertex[0].avertex().valsmooth<k1>(input_a, PT4_13_a_a2.vertex[0].tvertex());
    state_datatype PT4_K1p_0_nonladder = PT4_13_p_p2.vertex[0].pvertex().valsmooth<k1>(input_p, PT4_13_p_p2.vertex[0].pvertex());
    state_datatype PT4_K1t_0_nonladder_a = PT4_13_t_a2.vertex[0].tvertex().valsmooth<k1>(input_t, PT4_13_t_a2.vertex[0].avertex());
    state_datatype PT4_K1t_0_nonladder_t = PT4_13_t_t2.vertex[0].tvertex().valsmooth<k1>(input_t, PT4_13_t_t2.vertex[0].avertex());

    // K2 in PT3
    vec<state_datatype> PT3_K2a_0 (3);
    vec<state_datatype> PT3_K2p_0 (3);
    vec<state_datatype> PT3_K2t_0 (3);
    // K2 in PT4
    vec<state_datatype> PT4_K2a_0_p1 (3);
    vec<state_datatype> PT4_K2p_0_a1 (3);
    vec<state_datatype> PT4_K2a_0_t1 (3);
    vec<state_datatype> PT4_K2p_0_t1 (3);
    vec<state_datatype> PT4_K2a_0_a2 (3);
    vec<state_datatype> PT4_K2a_0_p2 (3);
    vec<state_datatype> PT4_K2a_0_t2 (3);
    vec<state_datatype> PT4_K2p_0_a2 (3);
    vec<state_datatype> PT4_K2p_0_p2 (3);
    vec<state_datatype> PT4_K2p_0_t2 (3);
    vec<state_datatype> PT4_K2t_0_a2 (3);
    vec<state_datatype> PT4_K2t_0_t2 (3);

#if MAX_DIAG_CLASS >= 2
#ifdef KELDYSH_FORMALISM
    for (int iK2=0; iK2<3; ++iK2) {
#else
      int iK2 = 0;
#endif
        input_a.iK = iK2s[iK2];
        input_p.iK = iK2s[iK2];
        input_t.iK = iK2s[iK2];
        // K2 in PT3
        PT3_K2a_0[iK2] = PT3_K2a.vertex[0].avertex().valsmooth<k2>(input_a, PT3_K2a.vertex[0].tvertex());
        PT3_K2p_0[iK2] = PT3_K2p.vertex[0].pvertex().valsmooth<k2>(input_p, PT3_K2p.vertex[0].pvertex());
        PT3_K2t_0[iK2] = PT3_K2t.vertex[0].tvertex().valsmooth<k2>(input_t, PT3_K2t.vertex[0].avertex());
        // K2 in PT4
        PT4_K2a_0_p1[iK2] = PT4_31_a_p1.vertex[0].avertex().valsmooth<k2>(input_a, PT4_31_a_p1.vertex[0].tvertex());
        PT4_K2p_0_a1[iK2] = PT4_31_p_a1.vertex[0].pvertex().valsmooth<k2>(input_p, PT4_31_p_a1.vertex[0].pvertex());
        PT4_K2a_0_t1[iK2] = PT4_31_a_t1.vertex[0].avertex().valsmooth<k2>(input_a, PT4_31_a_t1.vertex[0].tvertex());
        PT4_K2p_0_t1[iK2] = PT4_31_p_t1.vertex[0].pvertex().valsmooth<k2>(input_p, PT4_31_p_t1.vertex[0].pvertex());
        PT4_K2a_0_a2[iK2] = PT4_31_a_a2.vertex[0].avertex().valsmooth<k2>(input_a, PT4_31_a_a2.vertex[0].tvertex());
        PT4_K2a_0_p2[iK2] = PT4_31_a_p2.vertex[0].avertex().valsmooth<k2>(input_a, PT4_31_a_p2.vertex[0].tvertex());
        PT4_K2a_0_t2[iK2] = PT4_31_a_t2.vertex[0].avertex().valsmooth<k2>(input_a, PT4_31_a_t2.vertex[0].tvertex());
        PT4_K2p_0_a2[iK2] = PT4_31_p_a2.vertex[0].pvertex().valsmooth<k2>(input_p, PT4_31_p_a2.vertex[0].pvertex());
        PT4_K2p_0_p2[iK2] = PT4_31_p_p2.vertex[0].pvertex().valsmooth<k2>(input_p, PT4_31_p_p2.vertex[0].pvertex());
        PT4_K2p_0_t2[iK2] = PT4_31_p_t2.vertex[0].pvertex().valsmooth<k2>(input_p, PT4_31_p_t2.vertex[0].pvertex());
        PT4_K2t_0_a2[iK2] = PT4_31_t_a2.vertex[0].tvertex().valsmooth<k2>(input_t, PT4_31_t_a2.vertex[0].avertex());
        PT4_K2t_0_t2[iK2] = PT4_31_t_t2.vertex[0].tvertex().valsmooth<k2>(input_t, PT4_31_t_t2.vertex[0].avertex());
#ifdef KELDYSH_FORMALISM
    }
#endif

#ifdef KELDYSH_FORMALISM
    state_datatype PT3_K2_exact = -(1./2.) * (2. - M_PI*M_PI/4.) * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
#else
    state_datatype PT3_K2_exact = - (2. - M_PI*M_PI/4.) * glb_U * pow(glb_U / (M_PI * (glb_Gamma + Lambda) / 2.), 2);
    std::cout << "PT3 K2 exact: " << PT3_K2_exact << "\n";
    std::cout << "Computed value: " << PT3_K2t_0[iK2] << "\n";
#endif
#endif

    // K3 in PT4
#ifdef KELDYSH_FORMALISM
    input_a.iK = 5;
    input_p.iK = 5;
    input_t.iK = 5;
#endif
    state_datatype PT4_K3a_0;
    state_datatype PT4_K3p_0;
    state_datatype PT4_K3t_0_aa;
    state_datatype PT4_K3t_0_ap;
    state_datatype PT4_K3t_0_pa;

#if MAX_DIAG_CLASS == 3
    PT4_K3a_0 = PT4_22_a_pp.vertex[0].avertex().valsmooth<k3>(input_a, PT4_22_a_pp.vertex[0].tvertex());
    PT4_K3p_0 = PT4_22_p_aa.vertex[0].pvertex().valsmooth<k3>(input_p, PT4_22_p_aa.vertex[0].pvertex());
    PT4_K3t_0_aa = PT4_22_t_aa.vertex[0].tvertex().valsmooth<k3>(input_t, PT4_22_t_aa.vertex[0].avertex());
    PT4_K3t_0_ap = PT4_22_t_ap.vertex[0].tvertex().valsmooth<k3>(input_t, PT4_22_t_ap.vertex[0].avertex());
    PT4_K3t_0_pa = PT4_22_t_pa.vertex[0].tvertex().valsmooth<k3>(input_t, PT4_22_t_pa.vertex[0].avertex());
#endif

    // values to be printed to log
    vec<state_datatype> check_values {PT2_K1a_0 + PT2_K1p_0
                            , PT2_K1t_0
                            , (PT2_K1a_0 - PT2_K1_exact)/PT2_K1_exact
                            , (-PT2_K1p_0 - PT2_K1_exact)/PT2_K1_exact
                            , (PT3_K1a_0 - PT3_K1_exact)/PT3_K1_exact
                            , (PT3_K1p_0 - PT3_K1_exact)/PT3_K1_exact
                            , (PT3_K1t_0 - PT3_K1_exact)/PT3_K1_exact
#if MAX_DIAG_CLASS > 1
                            , (PT3_K2a_0[0] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[0] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[0] - PT3_K2_exact)/PT3_K2_exact
#ifdef KELDYSH_FORMALISM
                            , (PT3_K2a_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[1] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2a_0[2] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2p_0[2] - PT3_K2_exact)/PT3_K2_exact
                            , (PT3_K2t_0[2] - PT3_K2_exact)/PT3_K2_exact
#endif
                            , (PT4_K2a_0_p1[0] + PT4_K2p_0_a1[0])/(std::abs(PT4_K2a_0_p1[0]) + std::abs(PT4_K2p_0_a1[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2a_0_p1[1] + PT4_K2p_0_a1[1])/(std::abs(PT4_K2a_0_p1[1]) + std::abs(PT4_K2p_0_a1[1]))
                            , (PT4_K2a_0_p1[2] + PT4_K2p_0_a1[2])/(std::abs(PT4_K2a_0_p1[2]) + std::abs(PT4_K2p_0_a1[2]))
#endif
                            , (PT4_K2a_0_t1[0] + PT4_K2p_0_t1[0])/(std::abs(PT4_K2a_0_t1[0]) + std::abs(PT4_K2p_0_t1[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2a_0_t1[1] + PT4_K2p_0_t1[1])/(std::abs(PT4_K2a_0_t1[1]) + std::abs(PT4_K2p_0_t1[1]))
                            , (PT4_K2a_0_t1[2] + PT4_K2p_0_t1[2])/(std::abs(PT4_K2a_0_t1[2]) + std::abs(PT4_K2p_0_t1[2]))
#endif
                            , (PT4_K2a_0_a2[0] + PT4_K2p_0_p2[0])/(std::abs(PT4_K2a_0_a2[0]) + std::abs(PT4_K2p_0_p2[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2a_0_a2[1] + PT4_K2p_0_p2[1])/(std::abs(PT4_K2a_0_a2[1]) + std::abs(PT4_K2p_0_p2[1]))
                            , (PT4_K2a_0_a2[2] + PT4_K2p_0_p2[2])/(std::abs(PT4_K2a_0_a2[2]) + std::abs(PT4_K2p_0_p2[2]))
#endif
                            , (PT4_K2a_0_p2[0] + PT4_K2p_0_a2[0])/(std::abs(PT4_K2a_0_p2[0]) + std::abs(PT4_K2p_0_a2[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2a_0_p2[1] + PT4_K2p_0_a2[1])/(std::abs(PT4_K2a_0_p2[1]) + std::abs(PT4_K2p_0_a2[1]))
                            , (PT4_K2a_0_p2[2] + PT4_K2p_0_a2[2])/(std::abs(PT4_K2a_0_p2[2]) + std::abs(PT4_K2p_0_a2[2]))
#endif
                            , (PT4_K2a_0_t2[0] + PT4_K2p_0_t2[0])/(std::abs(PT4_K2a_0_t2[0]) + std::abs(PT4_K2p_0_t2[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2a_0_t2[1] + PT4_K2p_0_t2[1])/(std::abs(PT4_K2a_0_t2[1]) + std::abs(PT4_K2p_0_t2[1]))
                            , (PT4_K2a_0_t2[2] + PT4_K2p_0_t2[2])/(std::abs(PT4_K2a_0_t2[2]) + std::abs(PT4_K2p_0_t2[2]))
#endif
                            , (PT4_K2t_0_a2[0] + PT4_K2t_0_t2[0])/(std::abs(PT4_K2t_0_a2[0]) + std::abs(PT4_K2t_0_t2[0]))
#ifdef KELDYSH_FORMALISM
                            , (PT4_K2t_0_a2[1] + PT4_K2t_0_t2[1])/(std::abs(PT4_K2t_0_a2[1]) + std::abs(PT4_K2t_0_t2[1]))
                            , (PT4_K2t_0_a2[2] + PT4_K2t_0_t2[2])/(std::abs(PT4_K2t_0_a2[2]) + std::abs(PT4_K2t_0_t2[2]))
#endif
                            , (PT4_K3a_0 + PT4_K3p_0)/(std::abs(PT4_K3a_0) + std::abs(PT4_K3p_0))
                            , (PT4_K3t_0_aa + PT4_K3t_0_ap + PT4_K3t_0_pa)/(std::abs(PT4_K3t_0_aa) + std::abs(PT4_K3t_0_ap) + std::abs(PT4_K3t_0_pa))
#endif
                            };

    // print to log
    for (int i=0; i<check_values.size(); ++i) {
        print(check_labels[i], check_values[i], true);
    }

    print("----------------------", true);

    /*
    // Compute K1a contributions in PT4, using
    // (22):   K1a in PT2
    // (13_1): K1a in PT3
    // (13_2): K2a in PT3
    // (31_2): K2a in PT3
    State<Q> PT4_K1a22 (Lambda);
    State<Q> PT4_K1a13_1 (Lambda);
    State<Q> PT4_K1a13_2 (Lambda);
    State<Q> PT4_K1a31_2 (Lambda);

    t0 = get_time();
    bubble_function(PT4_K1a22.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false);
    bubble_function(PT4_K1a13_1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false);
    bubble_function(PT4_K1a13_2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false);
    bubble_function(PT4_K1a31_2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false);
    print("Computed K1 in PT4.", true);
    get_time(t0);

    // FDT checks
    print("Check K1a in PT4 (22):", true);
    check_FDTs(PT4_K1a22);
    print("Check K1a in PT4 (13_1):", true);
    check_FDTs(PT4_K1a13_1);
    print("Check K1a in PT4 (13_2):", true);
    check_FDTs(PT4_K1a13_2);
    print("Check K1a in PT4 (31_2):", true);
    check_FDTs(PT4_K1a31_2);
    // */
}

#if MAX_DIAG_CLASS == 3
/**
 * Test K3 dynamics by computing SE diagrams in PT4 using different PT4 vertices, which should all give the same result.
 */
template <typename Q>
void test_K3_dynamics_SE_PT4(double Lambda) {
    // Initialize a bare state
    State<Q> bare (Lambda);
    bare.initialize();

    // Initialize a bare propagator
    Propagator<Q> G(Lambda, bare.selfenergy, 'g');

    // Compute K1a and K1p in PT2
    State<Q> PT2_K1a (Lambda);
    State<Q> PT2_K1p (Lambda);

    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);

    // Compute K1a and K1p in PT3
    State<Q> PT3_K1a (Lambda);
    State<Q> PT3_K1p (Lambda);

    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT3_K1p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 'p', false);

    // Compute K3a, K1p (ladder), K3p, K1a (ladder) in PT4
    State<Q> PT4_22_a_pp (Lambda);
    State<Q> PT4_13_p_p1 (Lambda);
    State<Q> PT4_22_p_aa (Lambda);
    State<Q> PT4_13_a_a1 (Lambda);

    bubble_function(PT4_22_a_pp.vertex, PT2_K1p.vertex, PT2_K1p.vertex, G, G, 'a', false);
    bubble_function(PT4_13_p_p1.vertex, bare.vertex, PT3_K1p.vertex, G, G, 'p', false);
    bubble_function(PT4_22_p_aa.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'p', false);
    bubble_function(PT4_13_a_a1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false);

    // a-channel:
    // close K3a (single diagram: PT2_K1p - a-bubble - PT2_K1p)
    State<Q> SE_K3a (Lambda);
    loop(SE_K3a.selfenergy, PT4_22_a_pp.vertex, G, false);

    // close K1p ladder (use 1-3 vertex in p-channel, since it contains only p-ladder)
    State<Q> SE_K1p_ladder (Lambda);
    loop(SE_K1p_ladder.selfenergy, PT4_13_p_p1.vertex, G, false);

    // p-channel:
    // close K3p (single diagram: PT2_K1a - p-bubble - PT2_K1a)
    State<Q> SE_K3p (Lambda);
    loop(SE_K3p.selfenergy, PT4_22_p_aa.vertex, G, false);

    // close K1a ladder (use 1-3 vertex in a-channel, since it contains only a-ladder)
    State<Q> SE_K1a_ladder (Lambda);
    loop(SE_K1a_ladder.selfenergy, PT4_13_a_a1.vertex, G, false);

    write_hdf("SE_K3a_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K3a);
    write_hdf("SE_K1p_ladder_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K1p_ladder);
    write_hdf("SE_K3p_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K3p);
    write_hdf("SE_K1a_ladder_U" + std::to_string(1./((glb_Gamma+Lambda)/2.)) + ".h5", Lambda, 1, SE_K1a_ladder);
}
#endif



#if MAX_DIAG_CLASS >= 2
/**
 * Function to test correctness of K2a when calculating a susceptibility (a K1a-object) //Notice that the same calculation
 * can be performed in the p-channel.
 * @param Lambda : Scale at which the calculation is done.
 */
template <typename Q>
void test_K2_correctness(double Lambda){

    bool write_flag = true; //Write out results in a HDF5 file

    State<Q> bare (Lambda);   //Bare state
    bare.initialize();  //Initialize bare state

    Propagator<Q> G(Lambda, bare.selfenergy, 'g'); //Bare propagator
    Propagator<Q> S(Lambda, bare.selfenergy, 's'); //Bare single-scale propagator

    //Create states for K1-calculations
    State<Q> PT2_K1a (Lambda);
    State<Q> PT2_K1p (Lambda);
    State<Q> PT2_K1t (Lambda);

    //Save K1-bubbles in separate objects - SOPT
    double t0 = get_time();
    bubble_function(PT2_K1a.vertex, bare.vertex, bare.vertex, G, G, 'a', false);
    bubble_function(PT2_K1p.vertex, bare.vertex, bare.vertex, G, G, 'p', false);
    bubble_function(PT2_K1t.vertex, bare.vertex, bare.vertex, G, G, 't', false);
    get_time(t0);

    State<Q> PT2_SE_a (Lambda);
    State<Q> PT2_SE_p (Lambda);
    State<Q> PT2_SE_t (Lambda);
    State<Q> PT2_SE_p_1 (Lambda);
    State<Q> PT2_SE_p_4 (Lambda);
    State<Q> PT2_SE_p_5 (Lambda);

    loop(PT2_SE_a.selfenergy, PT2_K1a.vertex, S, true);
    loop(PT2_SE_p.selfenergy, PT2_K1p.vertex, S, true);
    loop(PT2_SE_t.selfenergy, PT2_K1t.vertex, S, true);
#ifdef DEBUG_MODE
    loop(PT2_SE_p_1.selfenergy, PT2_K1p.vertex, S, true, 1);
    loop(PT2_SE_p_4.selfenergy, PT2_K1p.vertex, S, true, 4);
    loop(PT2_SE_p_5.selfenergy, PT2_K1p.vertex, S, true, 5);
#endif

    State<Q> PT3_K2a (Lambda);    //Create state for K2a calculation
    State<Q> PT3_K2a_ia (Lambda);
    State<Q> PT3_K2a_ib (Lambda);
    State<Q> PT3_K2a_iia (Lambda);
    State<Q> PT3_K2a_iib (Lambda);
    State<Q> PT3_K2a_iva (Lambda);
    State<Q> PT3_K2a_ivb (Lambda);
    State<Q> PT3_K2a_t (Lambda);

    State<Q> PT3_K2p (Lambda);
    State<Q> PT3_K2t (Lambda);
    State<Q> PT3_K2t_a (Lambda);
    State<Q> PT3_K2t_p (Lambda);

    //Do appropriate calculation for K2a with K1p and K1t being fed back into the left vertex. Notice part = 'L' to ensure
    //that the correct contributions are added on both sides. - TOPT
    t0 = get_time();
    bubble_function(PT3_K2a.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false);   // K2a in PT3

#ifdef DEBUG_MODE
    bubble_function(PT3_K2a_ia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 9);
    bubble_function(PT3_K2a_iia.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 11, 14);
    bubble_function(PT3_K2a_iib.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 7, 9);
    bubble_function(PT3_K2a_iva.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 9, 6);
    bubble_function(PT3_K2a_ivb.vertex, PT2_K1p.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'a', false, 'L', 16, 16, 6, 15);
#endif

    bubble_function(PT3_K2a_t.vertex, PT2_K1t.vertex, bare.vertex, G, G, 'a', false);   // K2a in PT3
    //PT3_K2a = read_hdf("PT4_check_of_K2a_K2_switchedcc_adap_m3m9_g501_101_nI1501_state_PT3_K2a", 0, 1);

    bubble_function(PT3_K2p.vertex, PT2_K1a.vertex + PT2_K1t.vertex, bare.vertex, G, G, 'p', false);    // K2p  in PT3
    //bubble_function(PT3_K2t.vertex, PT2_K1a.vertex + PT2_K1p.vertex, bare.vertex, G, G, 't', false);    // K2t  in PT3
    bubble_function(PT3_K2t_a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 't', false);    // K2t  in PT3
    bubble_function(PT3_K2t_p.vertex, PT2_K1p.vertex, bare.vertex, G, G, 't', false);    // K2t  in PT3
    PT3_K2t.vertex = PT3_K2t_a.vertex + PT3_K2t_p.vertex;
    get_time(t0);

    // full K2 in PT3
    State<Q> PT3_K2 (Lambda);
    PT3_K2.vertex[0].avertex() = PT3_K2a_t.vertex[0].avertex();
    //PT3_K2.vertex[0].pvertex() = PT3_K2p.vertex[0].pvertex();
    PT3_K2.vertex[0].tvertex() = PT3_K2t_a.vertex[0].tvertex();

    State<Q> PT3_K2at (Lambda);

    // K2 contribution to self-energy flow
    State<Q> PT3_SE (Lambda);
    State<Q> PT3_SE_a (Lambda);
    State<Q> PT3_SE_p (Lambda);
    State<Q> PT3_SE_t (Lambda);
    State<Q> PT3_SE_t_1 (Lambda);
    State<Q> PT3_SE_t_4 (Lambda);
    State<Q> PT3_SE_t_5 (Lambda);
    State<Q> PT3_SE_t_a (Lambda);
    State<Q> PT3_SE_t_a_1 (Lambda);
    State<Q> PT3_SE_t_a_4 (Lambda);
    State<Q> PT3_SE_t_a_5 (Lambda);
    State<Q> PT3_SE_t_p (Lambda);
    State<Q> PT3_SE_t_p_1 (Lambda);
    State<Q> PT3_SE_t_p_4 (Lambda);
    State<Q> PT3_SE_t_p_5 (Lambda);

    loop(PT3_SE.selfenergy, PT3_K2.vertex, S, true);
    loop(PT3_SE_a.selfenergy, PT3_K2a.vertex, S, true);
    loop(PT3_SE_p.selfenergy, PT3_K2p.vertex, S, true);
    loop(PT3_SE_t.selfenergy, PT3_K2t.vertex, S, true);

    loop(PT3_SE_t_a.selfenergy, PT3_K2t_a.vertex, S, true);
    loop(PT3_SE_t_p.selfenergy, PT3_K2t_p.vertex, S, true);
#ifdef DEBUG_MODE
    loop(PT3_SE_t_1.selfenergy, PT3_K2t.vertex, S, true, 1);
    loop(PT3_SE_t_4.selfenergy, PT3_K2t.vertex, S, true, 4);
    loop(PT3_SE_t_5.selfenergy, PT3_K2t.vertex, S, true, 5);
    loop(PT3_SE_t_a_1.selfenergy, PT3_K2t_a.vertex, S, true, 1);
    loop(PT3_SE_t_a_4.selfenergy, PT3_K2t_a.vertex, S, true, 4);
    loop(PT3_SE_t_a_5.selfenergy, PT3_K2t_a.vertex, S, true, 5);
    loop(PT3_SE_t_p_1.selfenergy, PT3_K2t_p.vertex, S, true, 1);
    loop(PT3_SE_t_p_4.selfenergy, PT3_K2t_p.vertex, S, true, 4);
    loop(PT3_SE_t_p_5.selfenergy, PT3_K2t_p.vertex, S, true, 5);
#endif

    State<Q> PT3_K1a (Lambda);    //Create state to compare with K1a
    t0 = get_time();
    bubble_function(PT3_K1a.vertex, PT2_K1a.vertex, bare.vertex, G, G, 'a', false);
    get_time(t0);

    State<Q> PT123_a = bare + PT2_K1a + PT3_K1a + PT3_K2a;  //Create vertex of the right side of BSE

    State<Q> PT4_K1a22 (Lambda);
    State<Q> PT4_K1a13_1 (Lambda);
    State<Q> PT4_K1a13_2 (Lambda);
    State<Q> PT4_K1a31_2 (Lambda);
    State<Q> PT4_K1a13_2_11e (Lambda); // A
    State<Q> PT4_K1a13_2_21e (Lambda); // B
    State<Q> PT4_K1a13_2_11o (Lambda); // C
    State<Q> PT4_K1a13_2_21o (Lambda); // D
    State<Q> PT4_K1a13_2_12o (Lambda); // TST3TC D
    State<Q> PT4_K1a13_2_22o (Lambda); // F
    State<Q> PT4_K1a13_2_ia (Lambda);
    State<Q> PT4_K1a13_2_ib (Lambda);
    State<Q> PT4_K1a13_2_iia (Lambda);
    State<Q> PT4_K1a13_2_iib (Lambda);
    State<Q> PT4_K1a13_2_iva (Lambda);
    State<Q> PT4_K1a13_2_ivb (Lambda);

    t0 = get_time();
    bubble_function(PT4_K1a22.vertex, PT2_K1a.vertex, PT2_K1a.vertex, G, G, 'a', false);
    bubble_function(PT4_K1a13_1.vertex, bare.vertex, PT3_K1a.vertex, G, G, 'a', false);

    bubble_function(PT4_K1a13_2.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false);
    bubble_function(PT4_K1a31_2.vertex, PT3_K2a.vertex, bare.vertex, G, G, 'a', false);
    get_time(t0);
    t0 = get_time();
#ifdef DEBUG_MODE
    bubble_function(PT4_K1a13_2_11e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 0, 16, 16, 16); // A
    bubble_function(PT4_K1a13_2_21e.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 1, 16, 16, 16); // B
    bubble_function(PT4_K1a13_2_11o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 2, 16, 16, 16); // C
    bubble_function(PT4_K1a13_2_21o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 3, 16, 16, 16); // D
    bubble_function(PT4_K1a13_2_12o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 10, 16, 16, 16); // TST3TC D
    bubble_function(PT4_K1a13_2_22o.vertex, bare.vertex, PT3_K2a.vertex, G, G, 'a', false, 11, 16, 16, 16); // F

    bubble_function(PT4_K1a13_2_ia.vertex, bare.vertex, PT3_K2a_ia.vertex, G, G, 'a', false, 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_ib.vertex, bare.vertex, PT3_K2a_ib.vertex, G, G, 'a', false, 2, 9, 16, 16);
    bubble_function(PT4_K1a13_2_iia.vertex, bare.vertex, PT3_K2a_iia.vertex, G, G, 'a', false, 11, 6, 16, 16);
    bubble_function(PT4_K1a13_2_iib.vertex, bare.vertex, PT3_K2a_iib.vertex, G, G, 'a', false, 10, 11, 16, 16);
    bubble_function(PT4_K1a13_2_iva.vertex, bare.vertex, PT3_K2a_iva.vertex, G, G, 'a', false, 11, 15, 16, 16);
    bubble_function(PT4_K1a13_2_ivb.vertex, bare.vertex, PT3_K2a_ivb.vertex, G, G, 'a', false, 2, 9, 16, 16);
#endif
    get_time(t0);

    cvec K1a_diff(nBOS);
    for(int iw=0; iw<nBOS; ++iw){
        K1a_diff[iw] = PT4_K1a22.vertex[0].avertex().K1_val(0, iw, 0) - PT2_K1a.vertex[0].avertex().K1_val(0, iw, 0);
    }

    print("Testing correctness of K2a. Using U=" +std::to_string(glb_U)+ " and Lambda="+std::to_string(Lambda)
        +", the maximal difference between direct K1a and K1a over integration of K2a is " +std::to_string(K1a_diff.max_norm())+"." , true);
    if(write_flag) write_h5_rvecs("../Data/PT4_check_of_K2a_cleanup_GL_gW20_51_21_nI1501_U1", {"w",
                                                       "PT2_K1a_R", "PT2_K1a_I",
                                                       "PT2_K1p_R", "PT2_K1p_I",
                                                       "PT2_K1t_R", "PT2_K1t_I",
                                                       "PT2_SE_a_R", "PT2_SE_a_I",
                                                       "PT2_SE_p_R", "PT2_SE_p_I",
                                                       "PT2_SE_t_R", "PT2_SE_t_I",
                                                       "PT2_SE_p_1_R", "PT2_SE_p_1_I",
                                                       "PT2_SE_p_4_R", "PT2_SE_p_4_I",
                                                       "PT2_SE_p_5_R", "PT2_SE_p_5_I",
                                                       "PT3_K1a_R", "PT3_K1a_I",
                                                       "PT3_K2a_R", "PT3_K2a_I",
                                                       "PT3_K2a_ia_R", "PT3_K2a_ia_I",
                                                       "PT3_K2a_ib_R", "PT3_K2a_ib_I",
                                                       "PT3_K2a_iia_R", "PT3_K2a_iia_I",
                                                       "PT3_K2a_iib_R", "PT3_K2a_iib_I",
                                                       "PT3_K2a_iva_R", "PT3_K2a_iva_I",
                                                       "PT3_K2a_ivb_R", "PT3_K2a_ivb_I",
                                                       "PT3_K2a_t_R", "PT3_K2a_t_I",
                                                       "PT3_K2p_R", "PT3_K2p_I",
                                                       "PT3_K2t_R", "PT3_K2t_I",
                                                       "PT3_K2t_a_a_R", "PT3_K2t_a_a_I",
                                                       "PT3_K2t_a_p_R", "PT3_K2t_a_p_I",
                                                       "PT3_K2t_a_t_R", "PT3_K2t_a_t_I",
                                                       "PT3_K2t_p_a_R", "PT3_K2t_p_a_I",
                                                       "PT3_K2t_p_p_R", "PT3_K2t_p_p_I",
                                                       "PT3_K2t_p_t_R", "PT3_K2t_p_t_I",
                                                       "PT3_SE_R", "PT3_SE_I",
                                                       "PT3_SE_a_R", "PT3_SE_a_I",
                                                       "PT3_SE_p_R", "PT3_SE_p_I",
                                                       "PT3_SE_t_R", "PT3_SE_t_I",
                                                       "PT3_SE_t_1_R", "PT3_SE_t_1_I",
                                                       "PT3_SE_t_4_R", "PT3_SE_t_4_I",
                                                       "PT3_SE_t_5_R", "PT3_SE_t_5_I",
                                                       "PT3_SE_t_a_R", "PT3_SE_t_a_I",
                                                       "PT3_SE_t_a_1_R", "PT3_SE_t_a_1_I",
                                                       "PT3_SE_t_a_4_R", "PT3_SE_t_a_4_I",
                                                       "PT3_SE_t_a_5_R", "PT3_SE_t_a_5_I",
                                                       "PT3_SE_t_p_R", "PT3_SE_t_p_I",
                                                       "PT3_SE_t_p_1_R", "PT3_SE_t_p_1_I",
                                                       "PT3_SE_t_p_4_R", "PT3_SE_t_p_4_I",
                                                       "PT3_SE_t_p_5_R", "PT3_SE_t_p_5_I",
                                                       "PT4_K1a22_R", "PT4_K1a22_I",
                                                       "PT4_K1a13_1_R", "PT4_K1a13_1_I",
                                                       "PT4_K1a13_2_R", "PT4_K1a13_2_I",
                                                       "PT4_K1a31_2_R", "PT4_K1a31_2_I",
                                                       "PT4_K1a13_2_11e_R", "PT4_K1a13_2_11e_I",
                                                       "PT4_K1a13_2_21e_R", "PT4_K1a13_2_21e_I",
                                                       "PT4_K1a13_2_11o_R", "PT4_K1a13_2_11o_I",
                                                       "PT4_K1a13_2_21o_R", "PT4_K1a13_2_21o_I",
                                                       "PT4_K1a13_2_12o_R", "PT4_K1a13_2_12o_I",
                                                       "PT4_K1a13_2_22o_R", "PT4_K1a13_2_22o_I",
                                                       "PT4_K1a13_2_ia_R", "PT4_K1a13_2_ia_I",
                                                       "PT4_K1a13_2_ib_R", "PT4_K1a13_2_ib_I",
                                                       "PT4_K1a13_2_iia_R", "PT4_K1a13_2_iia_I",
                                                       "PT4_K1a13_2_iib_R", "PT4_K1a13_2_iib_I",
                                                       "PT4_K1a13_2_iva_R", "PT4_K1a13_2_iva_I",
                                                       "PT4_K1a13_2_ivb_R", "PT4_K1a13_2_ivb_I"
                                                       },
                                  {bfreqs,
                                   PT2_K1a.vertex[0].avertex().K1.real(), PT2_K1a.vertex[0].avertex().K1.imag(),
                                   PT2_K1p.vertex[0].pvertex().K1.real(), PT2_K1p.vertex[0].pvertex().K1.imag(),
                                   PT2_K1t.vertex[0].tvertex().K1.real(), PT2_K1t.vertex[0].tvertex().K1.imag(),
                                   PT2_SE_a.selfenergy.Sigma.real(), PT2_SE_a.selfenergy.Sigma.imag(),
                                   PT2_SE_p.selfenergy.Sigma.real(), PT2_SE_p.selfenergy.Sigma.imag(),
                                   PT2_SE_t.selfenergy.Sigma.real(), PT2_SE_t.selfenergy.Sigma.imag(),
                                   PT2_SE_p_1.selfenergy.Sigma.real(), PT2_SE_p_1.selfenergy.Sigma.imag(),
                                   PT2_SE_p_4.selfenergy.Sigma.real(), PT2_SE_p_4.selfenergy.Sigma.imag(),
                                   PT2_SE_p_5.selfenergy.Sigma.real(), PT2_SE_p_5.selfenergy.Sigma.imag(),
                                   PT3_K1a.vertex[0].avertex().K1.real(), PT3_K1a.vertex[0].avertex().K1.imag(),
                                   PT3_K2a.vertex[0].avertex().K2.real(), PT3_K2a.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_ia.vertex[0].avertex().K2.real(), PT3_K2a_ia.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_ib.vertex[0].avertex().K2.real(), PT3_K2a_ib.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_iia.vertex[0].avertex().K2.real(), PT3_K2a_iia.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_iib.vertex[0].avertex().K2.real(), PT3_K2a_iib.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_iva.vertex[0].avertex().K2.real(), PT3_K2a_iva.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_ivb.vertex[0].avertex().K2.real(), PT3_K2a_ivb.vertex[0].avertex().K2.imag(),
                                   PT3_K2a_t.vertex[0].avertex().K2.real(), PT3_K2a_t.vertex[0].avertex().K2.imag(),
                                   PT3_K2p.vertex[0].avertex().K2.real(), PT3_K2p.vertex[0].avertex().K2.imag(),
                                   PT3_K2t.vertex[0].avertex().K2.real(), PT3_K2t.vertex[0].avertex().K2.imag(),
                                   PT3_K2t_a.vertex[0].avertex().K2.real(), PT3_K2t_a.vertex[0].avertex().K2.imag(),
                                   PT3_K2t_a.vertex[0].pvertex().K2.real(), PT3_K2t_a.vertex[0].pvertex().K2.imag(),
                                   PT3_K2t_a.vertex[0].tvertex().K2.real(), PT3_K2t_a.vertex[0].tvertex().K2.imag(),
                                   PT3_SE.selfenergy.Sigma.real(), PT3_SE.selfenergy.Sigma.imag(),
                                   PT3_SE_a.selfenergy.Sigma.real(), PT3_SE_a.selfenergy.Sigma.imag(),
                                   PT3_SE_p.selfenergy.Sigma.real(), PT3_SE_p.selfenergy.Sigma.imag(),
                                   PT3_SE_t.selfenergy.Sigma.real(), PT3_SE_t.selfenergy.Sigma.imag(),
                                   PT3_SE_t_1.selfenergy.Sigma.real(), PT3_SE_t_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_4.selfenergy.Sigma.real(), PT3_SE_t_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_5.selfenergy.Sigma.real(), PT3_SE_t_5.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a.selfenergy.Sigma.real(), PT3_SE_t_a.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_1.selfenergy.Sigma.real(), PT3_SE_t_a_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_4.selfenergy.Sigma.real(), PT3_SE_t_a_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_a_5.selfenergy.Sigma.real(), PT3_SE_t_a_5.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p.selfenergy.Sigma.real(), PT3_SE_t_p.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_1.selfenergy.Sigma.real(), PT3_SE_t_p_1.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_4.selfenergy.Sigma.real(), PT3_SE_t_p_4.selfenergy.Sigma.imag(),
                                   PT3_SE_t_p_5.selfenergy.Sigma.real(), PT3_SE_t_p_5.selfenergy.Sigma.imag(),
                                   PT4_K1a22.vertex[0].avertex().K1.real(), PT4_K1a22.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_1.vertex[0].avertex().K1.real(), PT4_K1a13_1.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2.vertex[0].avertex().K1.real(), PT4_K1a13_2.vertex[0].avertex().K1.imag(),
                                   PT4_K1a31_2.vertex[0].avertex().K1.real(), PT4_K1a31_2.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_11e.vertex[0].avertex().K1.real(), PT4_K1a13_2_11e.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_21e.vertex[0].avertex().K1.real(), PT4_K1a13_2_21e.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_11o.vertex[0].avertex().K1.real(), PT4_K1a13_2_11o.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_21o.vertex[0].avertex().K1.real(), PT4_K1a13_2_21o.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_12o.vertex[0].avertex().K1.real(), PT4_K1a13_2_12o.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_22o.vertex[0].avertex().K1.real(), PT4_K1a13_2_22o.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_ia.vertex[0].avertex().K1.real(), PT4_K1a13_2_ia.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_ib.vertex[0].avertex().K1.real(), PT4_K1a13_2_ib.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_iia.vertex[0].avertex().K1.real(), PT4_K1a13_2_iia.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_iib.vertex[0].avertex().K1.real(), PT4_K1a13_2_iib.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_iva.vertex[0].avertex().K1.real(), PT4_K1a13_2_iva.vertex[0].avertex().K1.imag(),
                                   PT4_K1a13_2_ivb.vertex[0].avertex().K1.real(), PT4_K1a13_2_ivb.vertex[0].avertex().K1.imag()});

    //write_hdf("PT4_check_of_K2a_K2_switchedcc_t_update_symmrev_new11_SE_symm_full_adap_m3m9_gW10_501_101_nI1501_U1_state_PT3_K2a", 0, 1, PT3_K2a);
}

/**
 * Master function to test both consistency and correctness of K2-class
 * @param Lambda
 */
template<typename Q>
void test_K2(double Lambda, bool test_consistency){


    //First test consistency
    if(test_consistency) {
        bool K2a = test_K2_consistency<Q>(Lambda, 'a');    //Consistency of a-channel
        bool K2p = test_K2_consistency<Q>(Lambda, 'p');    //Consistency of p-channel
        bool K2t = test_K2_consistency<Q>(Lambda, 't');    //Consistency of t-channel

        if(K2a&&K2p&&K2t)
            test_K2_correctness<Q>(Lambda);
    }

    test_K2_correctness<Q>(Lambda);

}

/**
 * test-integrand for below function test_integrate_over_K1()
 */
template <typename Q>
class TestIntegrandK1a{
    public:
        const double Lambda;
        const double w, v, vp;
        const char channel;
        State<Q> SOPTstate;
        TestIntegrandK1a(double Lambda_in, char channel_in, double w, double v, double vp)
        : Lambda(Lambda_in), channel(channel_in), w(w), v(v), vp(vp) {
            SOPTstate = State<Q>(Lambda);
            SOPTstate.initialize();             // initialize state

            sopt_state(SOPTstate, Lambda);
        }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    void save_state() {
            write_hdf("SOPT_state.h5", 1.8, 1, SOPTstate);
        }

    auto operator() (double vpp) const -> Q {
            VertexInput input(0, vpp, v, vp + vpp, 0, 0, channel);
            return SOPTstate.vertex[0].half1().avertex.value(input, this->SOPTstate.vertex[0].half1().avertex) ;
            //return vpp*vpp;
        }
    };


/**
 * test function; can be used to test the integration and the interpolation
 * @tparam Q
 * @param Lambda
 */
template <typename Q>
void test_integrate_over_K1(double Lambda) {
    double v = 1e2;
    double vp = -1e2;
    TestIntegrandK1a<Q> IntegrandK1a(Lambda, 'a', 0., v, vp);
    IntegrandK1a.save_integrand();
    double exact_result = -glb_U*glb_U/4;

    double vmax = 1e2;
    //Integrand& integrand, const double vmin, const double vmax, double w_half, const vec<double>& freqs, const double Delta, const int num_freqs;
    vec<double> freqs = {};
    int num_freqs = 0;
    double Delta = (glb_Gamma+Lambda)/2.;
    Q result = integrator<Q,0>(IntegrandK1a, -vmax, vmax, 0, freqs, Delta)/ (2*M_PI);
    print("Result of the integration over K1a:", result, true);
    print("relative difference: ", (result-exact_result)/exact_result, true);


    TestIntegrandK1a<Q> IntegrandK1a2(Lambda, 't', (v-vp)*2, v, vp);
    Q result2 = integrator<Q,0>(IntegrandK1a2, -vmax, vmax, (v-vp), freqs, Delta)/ (2*M_PI);
    print("Result of the integration over K1a from the t-channel:", result2, true);
    print("relative difference: ", (result2-exact_result)/exact_result, true);


    TestIntegrandK1a<Q> IntegrandK1a3(Lambda, 'p', (v+vp)*2, v, vp);
    Q result3 = integrator<Q,0>(IntegrandK1a3, -vmax, vmax, (v+vp), freqs, Delta)/ (2*M_PI);
    print("Result of the integration over K1a from p-channel:", result3, true);
    print("relative difference: ", (result3-exact_result)/exact_result, true);


    IntegrandK1a2.save_integrand(vmax);
    IntegrandK1a3.save_integrand(vmax);
    IntegrandK1a.save_state();
}

#if not defined(KELDYSH_FORMALISM) and defined(ZERO_TEMP)
auto SOPT_K1a(double w, double Lambda) -> double {
    double Delta = (glb_Gamma + Lambda) / 2.;
    double result;
    if (w == 0.)
        result = - glb_U*glb_U * Delta / M_PI / (Delta*Delta + glb_mu*glb_mu);
    else
        result = - glb_U*glb_U * Delta/M_PI / (std::std::abs(w)*(2*Delta + std::std::abs(w))) * log(1 + (std::std::abs(w)*(2*Delta + std::std::abs(w)))/(Delta*Delta + glb_mu*glb_mu));
    return result;
}
auto SOPT_K1a_diff(double w, double Lambda) -> double {
    double Delta = (glb_Gamma + Lambda) / 2.;
    if (w == 0.) return - glb_U*glb_U      /2 / M_PI / (Delta*Delta + glb_mu*glb_mu)
                        + glb_U*glb_U * Delta / M_PI / (Delta*Delta + glb_mu*glb_mu) / (Delta*Delta + glb_mu*glb_mu) * Delta ;
    double term1 =  - glb_U*glb_U      / 2/M_PI / (std::std::abs(w)*(2*Delta + std::std::abs(w))) * log(1 + (std::std::abs(w)*(2*Delta + std::std::abs(w)))/(Delta*Delta + glb_mu*glb_mu));
    double term2 =  + glb_U*glb_U * Delta /M_PI / (std::std::abs(w)*(2*Delta + std::std::abs(w))) * log(1 + (std::std::abs(w)*(2*Delta + std::std::abs(w)))/(Delta*Delta + glb_mu*glb_mu)) / (2*Delta + std::std::abs(w));
    double term3 =  - glb_U*glb_U * Delta /M_PI / (std::std::abs(w)*(2*Delta + std::std::abs(w)))     /(   1 + (std::std::abs(w)*(2*Delta + std::std::abs(w)))/(Delta*Delta + glb_mu*glb_mu))
            *(
                    std::std::abs(w)                      / (Delta*Delta + glb_mu*glb_mu)
                   -std::std::abs(w) * (2*Delta + std::std::abs(w)) / (Delta*Delta + glb_mu*glb_mu) / (Delta*Delta + glb_mu*glb_mu) * Delta
                    );
    return term1 + term2 + term3;
}

/**
 * integrand for selfenergy in TOPT (SIAM)
 */
template <typename Q>
class Integrand_TOPT_SE{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    Propagator<Q>& barePropagator;
    Integrand_TOPT_SE(double Lambda_in, double w, double v, bool diff, Propagator<Q>& barePropagator)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), barePropagator(barePropagator) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(- v + vpp, Lambda) * barePropagator.valsmooth(0, vpp, 0) ;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class Integrand_TOPTK2a{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    Integrand_TOPTK2a(double Lambda_in, double w, double v, bool diff, Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * glb_U ;
        //return vpp*vpp;
    }
};

/**
 * integrand for K3a in FOPT (SIAM)
 */
template <typename Q>
class Integrand_FOPTK3a{
public:
    const double Lambda;
    const double w, v, vp;
    const bool diff;
    const char channel = 'a';
    Bubble<Q> Pi;
    Integrand_FOPTK3a(double Lambda_in, double w, double v, double vp, bool diff, Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) { }

    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * SOPT_K1a(vp + vpp, Lambda) ;
        //return vpp*vpp;
    }
};

template <typename Q>
void test_PT_state(std::string outputFileName, double Lambda, bool diff) {
    double vmax = 100;
    double Delta = (glb_Gamma + Lambda) / 2.;
    State<Q> bareState (Lambda);
    bareState.initialize();  //a state with a bare vertex and a self-energy initialized at the Hartree value
    Propagator<Q> barePropagator(Lambda, bareState.selfenergy, 'g');    //Bare propagator
    Bubble<Q> Pi(barePropagator, barePropagator, diff);

    State<state_datatype> PT_state(Lambda);

    for (int i = 0; i<nFER-0; i++) {
        double v = PT_state.selfenergy.frequencies.ws[i];
        Integrand_TOPT_SE<Q> IntegrandSE(Lambda, 0, v, diff, barePropagator);
        Q val_SE = 1./(2*M_PI) * integrator<Q,1>(IntegrandSE, -vmax, vmax, std::std::abs(0.), {v}, Delta, false);
        PT_state.selfenergy.setself(0, i, 0, val_SE);
    }

    for (int i = 0; i<nBOS-0; i++) {
        double w = PT_state.vertex[0].avertex().frequencies.b_K1.ws[i];
        Q val_K1;
        if (diff) val_K1 = SOPT_K1a_diff(w, Lambda);
        else val_K1 = SOPT_K1a(w, Lambda);
        PT_state.vertex[0].avertex().K1_setvert(0, i, 0, val_K1);
        PT_state.vertex[0].pvertex().K1_setvert(0, i, 0, -val_K1);
        PT_state.vertex[0].tvertex().K1_setvert(0, i, 0, -val_K1*val_K1/glb_U);
    }

#if MAX_DIAG_CLASS > 1
    for (int i = 0; i<nBOS2-0; i++) {
        for (int j = 1; j<nFER2-1; j++) {
            double w = PT_state.vertex[0].avertex().frequencies.b_K2.ws[i];
            double v = PT_state.vertex[0].avertex().frequencies.f_K2.ws[j];
            Integrand_TOPTK2a<Q> IntegrandK2(Lambda, w, v, diff, Pi);
            Q val_K2 = 1./(2*M_PI) * integrator<Q,3>(IntegrandK2, -vmax, vmax, std::std::abs(w/2), {v, w+v, w-v}, Delta, true);
            PT_state.vertex[0].avertex().K2_setvert(0, i, j, 0, val_K2);
            PT_state.vertex[0].pvertex().K2_setvert(0, i, j, 0, val_K2);
            PT_state.vertex[0].tvertex().K2_setvert(0, i, j, 0, val_K2);
        }
    }
#endif
#if MAX_DIAG_CLASS > 2
    for (int i = 0; i<nBOS3-0; i++) {
        for (int j = 1; j<nFER3-1; j++) {
            for (int k = 1; k<nFER3-1; k++) {
                double w = PT_state.vertex[0].avertex().frequencies.b_K3.ws[i];
                double v = PT_state.vertex[0].avertex().frequencies.f_K3.ws[j];
                double vp= PT_state.vertex[0].avertex().frequencies.f_K3.ws[k];
                Integrand_FOPTK3a<Q> IntegrandK3(Lambda, w, v, vp, diff, Pi);
                Q val_K3 = 1./(2*M_PI) * integrator<Q,6>(IntegrandK3, -vmax, vmax, std::std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                PT_state.vertex[0].avertex().K3_setvert(0, i, j, k, 0, val_K3);
                PT_state.vertex[0].pvertex().K3_setvert(0, i, j, k, 0, -val_K3);
                Integrand_FOPTK3a<Q> IntegrandK3_ap(Lambda, w, -v, vp, diff, Pi);
                Q val_K3_ap = 1./(2*M_PI) * integrator<Q,6>(IntegrandK3_ap, -vmax, vmax, std::std::abs(w/2), {v, vp, w+v, w-v, w+vp, w-vp}, Delta, true);
                PT_state.vertex[0].tvertex().K3_setvert(0, i, j, k, 0, -2*(val_K3-val_K3_ap));
            }
        }
    }
#endif

    write_hdf(outputFileName + "_exact", Lambda, 1, PT_state);

    State<Q> state_cpp (Lambda);   // create final and initial state
    state_cpp.initialize();             // initialize state

    fopt_state(state_cpp, Lambda);

    write_hdf(outputFileName + "_cpp", Lambda, 1, state_cpp);


    State<Q> state_diff = state_cpp - PT_state;

    write_hdf(outputFileName + "_diff", Lambda, 1, state_diff);
    print("K1-difference: ", state_diff.vertex[0].half1().norm_K1(0), true);
#if MAX_DIAG_CLASS > 1
    print("K2-difference: ", state_diff.vertex[0].half1().norm_K2(0), true);
#endif
#if MAX_DIAG_CLASS > 2
    print("K3-difference: ", state_diff.vertex[0].half1().norm_K3(0), true);
#endif

}





/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class K1rdot_PIa_K1p_exact_K2{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    const char channel = 'a';
    const bool diff;

    const Bubble<Q> Pi;
    K1rdot_PIa_K1p_exact_K2(double Lambda_in, double w, double v, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a_diff(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * glb_U ;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class K1rdot_PIa_K1p_exact_K3{
public:
    const double Lambda;
    const double w, v, vp;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    K1rdot_PIa_K1p_exact_K3(double Lambda_in, double w, double v, double vp, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) { }


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {
        return SOPT_K1a_diff(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * SOPT_K1a(vp + vpp, Lambda) ;
        //return vpp*vpp;
    }
};


/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K1{
public:
    const double Lambda;
    const double w;
    const double v = 0.;
    const double vp = 0.;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    const Bubble<Q> Pi;
    IntegranddGammaC_exact_K1(double Lambda_in, double w, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, vpp, false, Pi);
        state_datatype val_K2 = 1./(2*M_PI) * integrator<state_datatype,2>(IntegrandK2, -vmax, vmax, std::std::abs(w/2), {vpp, vp}, Delta);
        return -glb_U * Pi.value(0, w, vpp, 0, 'a') * val_K2;
        //return vpp*vpp;
    }
};

/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K2{
public:
    const double Lambda;
    const double w, v;
    const double vp = 0.;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    IntegranddGammaC_exact_K2(double Lambda_in, double w, double v, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, vpp, false, Pi);
        state_datatype val_K2 = 1./(2*M_PI) * integrator<state_datatype,1>(IntegrandK2, -vmax, vmax, std::std::abs(w/2), {vpp}, Delta);
        return -SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * val_K2;
        //return vpp*vpp;
    }
};


/**
 * integrand for K2a in TOPT (SIAM)
 */
template <typename Q>
class IntegranddGammaC_exact_K3{
public:
    const double Lambda;
    const double w, v, vp;
    double vmax = 1e3;
    const double Delta = (Lambda * glb_Gamma)/2;
    const char channel = 'a';
    const bool diff;

    Bubble<Q> Pi;
    IntegranddGammaC_exact_K3(double Lambda_in, double w, double v, double vp, bool diff, const Bubble<Q>& Pi)
            : Lambda(Lambda_in), w(w), v(v), vp(vp), diff(diff), Pi(Pi) {}


    void save_integrand(double vmax) {
        int npoints = 1e5;
        rvec freqs (npoints);

        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        double wl = -vmax;
        double wu = vmax;
        for (int i=0; i<npoints; ++i) {
            double vpp = wl + i * (wu-wl)/(npoints-1);
            Q integrand_value = (*this)(vpp);
            freqs[i] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i] = integrand_value;
            integrand_im[i] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }
    void save_integrand() {
        int npoints =nBOS*2-1;
        rvec freqs (npoints);
        double frac = 0.5;
        rvec integrand_re (npoints);
        rvec integrand_im (npoints);
        for (int i=0; i<nBOS-1; ++i) {
            double vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i];
            Q integrand_value = (*this)(vpp);
            freqs[i*2] = vpp;


#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2] = integrand_value;
            integrand_im[i*2] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif

            vpp = this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i] * (frac) + this->SOPTstate.vertex[0].half1().avertex.frequencies.b_K1.ws[i+1] * (1-frac);
            integrand_value = (*this)(vpp);
            freqs[i*2+1] = vpp;

#if defined(PARTICLE_HOLE_SYMM) and not defined(KELDYSH_FORMALISM)
            integrand_re[i*2+1] = integrand_value;
            integrand_im[i*2+1] = 0.;
#else
            integrand_re[i] = integrand_value.real();
            integrand_im[i] = integrand_value.imag();
#endif
        }

        std::string filename = "../Data/integrand_K1";
        filename += channel;
        filename += "_w=" + std::to_string(w) +"_v" + std::to_string(v) +  "_vp" + std::to_string(vp) + ".h5";
        write_h5_rvecs(filename,
                       {"v", "integrand_re", "integrand_im"},
                       {freqs, integrand_re, integrand_im});
    }

    auto operator() (double vpp) const -> Q {

        K1rdot_PIa_K1p_exact_K3<state_datatype> IntegrandK3(Lambda, w, vpp, vp, false, Pi);
        state_datatype val_K3 = 1./(2*M_PI) * integrator<state_datatype,3>(IntegrandK3, -vmax, vmax, std::std::abs(w/2), {vpp, vp, std::std::abs(vpp)-std::std::abs(vp)}, Delta);
        return -SOPT_K1a(v + vpp, Lambda) * Pi.value(0, w, vpp, 0, 'a') * val_K3;
        //return vpp*vpp;
    }
};



// to check central part of multi-loop flow equations:
// compute diagrams with non-symmetric intermediate results
void compute_non_symmetric_diags(const double Lambda, bool write_flag = false, int version=1, bool compute_exact=false) {
    State<state_datatype> bare (Lambda); // bare state
    bare.initialize();         // initialize bare state

    Propagator<state_datatype> G (Lambda, bare.selfenergy, 'g'); // bare propagator
    Propagator<state_datatype> S (Lambda, bare.selfenergy, 's'); // bare differentiated propagator = single scale propagator
    const Bubble<state_datatype> Pi(G, S, false);
    double Delta = (glb_Gamma + Lambda)/2;

    // Psi := K1p in PT2 + bare vertex
    State<state_datatype> Psi (Lambda);
    Psi.initialize();         // initialize bare state
    bubble_function(Psi.vertex, bare.vertex, bare.vertex, G, G, 'p', false);

    // K1a_dot in PT2
    State<state_datatype> PT2_K1adot (Lambda);
    bubble_function(PT2_K1adot.vertex, bare.vertex, bare.vertex, G, S, 'a', true);
    // K1p_dot in PT2
    State<state_datatype> PT2_K1pdot (Lambda);
    bubble_function(PT2_K1pdot.vertex, bare.vertex, bare.vertex, G, S, 'p', true);

    if (write_flag) {
        write_hdf("Psi_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, Psi);
        write_hdf("PT2_K1a_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1adot);
        write_hdf("PT2_K1p_dot_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, PT2_K1pdot);
    }

    std::vector<State<state_datatype>> central_bubblestates = {PT2_K1adot, PT2_K1pdot};

    //for (int i = 0; i < 2; i++){
    int i = version;
        State<state_datatype> centralstate_dot = central_bubblestates[i];

        // intermediate results
        State<state_datatype> K1rdot_PIa_K1p (Lambda);
        bubble_function(K1rdot_PIa_K1p.vertex, centralstate_dot.vertex, Psi.vertex, G, G, 'a', false);


        State<state_datatype> K1p_PIa_K1rdot (Lambda);
        bubble_function(K1p_PIa_K1rdot.vertex, Psi.vertex, centralstate_dot.vertex, G, G, 'a', false);


        if (write_flag) {
            write_hdf("K1rdot_PIa_K1p_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_hdf("K1p_PIa_K1rdot_UNREORDERED_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
        }

        Vertex<state_datatype> dGammaL_half1 = K1rdot_PIa_K1p.vertex;
        Vertex<state_datatype> dGammaR_half1 = K1p_PIa_K1rdot.vertex;
        dGammaL_half1[0].half1().reorder_due2antisymmetry(dGammaR_half1[0].half1());
        dGammaR_half1[0].half1().reorder_due2antisymmetry(dGammaL_half1[0].half1());
        K1rdot_PIa_K1p.vertex = dGammaL_half1;
        K1p_PIa_K1rdot.vertex = dGammaR_half1;

        // create non-symmetric vertex with differentiated vertex on the left
        GeneralVertex<state_datatype , non_symmetric> dGammaL(n_spin, Lambda);
        dGammaL[0].half1()  = dGammaL_half1[0].half1();  // assign half 1 to dGammaL
        dGammaL[0].half2() = dGammaR_half1[0].half1();  // assign half 2 as half 1 to dGammaR [symmetric -> left()=right()]

        // insert this non-symmetric vertex on the right of the bubble
        State<state_datatype> dGammaC_r(Lambda);
        bubble_function(dGammaC_r.vertex, Psi.vertex, dGammaL, G, G, 'a', false);


        // create non-symmetric vertex with differentiated vertex on the right (full dGammaR, containing half 1 and 2)
        GeneralVertex<state_datatype , non_symmetric> dGammaR (n_spin, Lambda);
        dGammaR[0].half1() = dGammaR_half1[0].half1();  // assign half 1
        dGammaR[0].half2() = dGammaL_half1[0].half1();  // assign half 2 as half 1 of dGammaL

        // insert this non-symmetric vertex on the left of the bubble
        State<state_datatype> dGammaC_l(Lambda);
        bubble_function(dGammaC_l.vertex, dGammaR, Psi.vertex, G, G, 'a', false);


        if (write_flag) {
            write_hdf("K1rdot_PIa_K1p_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1rdot_PIa_K1p);
            write_hdf("K1p_PIa_K1rdot_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, K1p_PIa_K1rdot);
            write_hdf("dGammaC_r_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_r);
            write_hdf("dGammaC_l_version" + std::to_string(i) + "_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5", Lambda, 1, dGammaC_l);
        }
    //}

    if (compute_exact) {
        /// Now computing vertex for version 1 (K1rdot = K1pdot):
        double vmax = 1e3;
        State<state_datatype> K1pdot_exact(Lambda);        // intermediate result: contains K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(K1pdot_exact, vmax, Delta)
        for (int iflat = 0; iflat < (nBOS - 2); iflat++) {
            int i = 1 + iflat;
            //for (int i = 1; i<nBOS2-1; i++) {
            //    for (int j = 1; j<nFER2-1; j++) {
            double w = K1pdot_exact.vertex[0].avertex().frequencies.b_K1.ws[i];
            state_datatype val_K1 = -SOPT_K1a_diff(w, Lambda);
            K1pdot_exact.vertex[0].pvertex().K1_setvert(0, i, 0, val_K1);
            //    }
        }
        write_hdf("K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  K1pdot_exact);

        State<state_datatype> K1rdot_diff = PT2_K1pdot - K1pdot_exact;        // intermediate result: contains K2 and K3
        write_hdf("K1rdot_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  K1rdot_diff);


        State<state_datatype> K1rdot_PIa_K1p_exact(Lambda);        // intermediate result: contains K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(K1rdot_PIa_K1p_exact, vmax, Delta)
        for (int iflat = 0; iflat < (nBOS2 - 2) * (nFER2 - 2); iflat++) {
            int i = 1 + iflat / (nFER2 - 2);
            int j = 1 + iflat - (i - 1) * (nFER2 - 2);
            //for (int i = 1; i<nBOS2-1; i++) {
            //    for (int j = 1; j<nFER2-1; j++) {
            double w = K1rdot_PIa_K1p_exact.vertex[0].avertex().frequencies.b_K2.ws[i];
            double v = K1rdot_PIa_K1p_exact.vertex[0].avertex().frequencies.f_K2.ws[j];
            K1rdot_PIa_K1p_exact_K2<state_datatype> IntegrandK2(Lambda, w, v, false, Pi);
            state_datatype val_K2 =
                    1. / (2 * M_PI) * integrator<state_datatype, 1>(IntegrandK2, -vmax, vmax, std::std::abs(w / 2), {v}, Delta);
            K1rdot_PIa_K1p_exact.vertex[0].avertex().K2_setvert(0, i, j, 0, val_K2);
            //    }
        }

#if MAX_DIAG_CLASS>2
#pragma omp parallel for schedule(dynamic) default(none) shared(K1rdot_PIa_K1p_exact, vmax, Delta)
        for (int iflat = 0; iflat < (nBOS3 - 2) * (nFER3 - 2) * (nFER3 - 2); iflat++) {
            int i = 1 + iflat / (nFER3 - 2) / (nFER3 - 2);
            int j = 1 + iflat / (nFER3 - 2) - (i - 1) * (nFER3 - 2);
            int k = 1 + iflat - (i - 1) * (nFER3 - 2) * (nFER3 - 2) - (j - 1) * (nFER3 - 2);
            //for (int i = 1; i<nBOS3-1; i++) {
            //    for (int j = 1; j<nFER3-1; j++) {
            //        for (int k = 1; k<nFER3-1; k++) {
            double w = K1rdot_PIa_K1p_exact.vertex[0].avertex().frequencies.b_K3.ws[i];
            double v = K1rdot_PIa_K1p_exact.vertex[0].avertex().frequencies.f_K3.ws[j];
            double vp = K1rdot_PIa_K1p_exact.vertex[0].avertex().frequencies.f_K3.ws[k];
            K1rdot_PIa_K1p_exact_K3<state_datatype> IntegrandK3(Lambda, w, v, vp, false, Pi);
            state_datatype val_K3 = 1. / (2 * M_PI) *
                                    integrator<state_datatype, 6>(IntegrandK3, -vmax, vmax, std::std::abs(w / 2),
                                                                  {v, vp, std::std::abs(w) - std::std::abs(vp), std::std::abs(w) + std::std::abs(vp),
                                                                   std::std::abs(w) - std::std::abs(v), std::std::abs(w) + std::std::abs(v)}, Delta);
            K1rdot_PIa_K1p_exact.vertex[0].avertex().K3_setvert(0, i, j, k, 0, val_K3);
            K1rdot_PIa_K1p_exact.vertex[0].pvertex().K3_setvert(0, i, j, k, 0, val_K3);
            K1rdot_PIa_K1p_exact.vertex[0].tvertex().K3_setvert(0, i, j, k, 0, val_K3);
            //        }
            //    }
        }
#endif
        write_hdf("K1rdot_PIa_K1p_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  K1rdot_PIa_K1p_exact);

        State<state_datatype> K1rdot_PIa_K1p_diff =
                K1rdot_PIa_K1p - K1rdot_PIa_K1p_exact;        // intermediate result: contains K2 and K3
        write_hdf("K1rdot_PIa_K1p_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  K1rdot_PIa_K1p_diff);


        State<state_datatype> dGammaC_exact(Lambda);        // final state: contains K1, K2 and K3

#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta)
        for (int i = 1; i < nBOS - 1; i++) {
            double w = dGammaC_exact.vertex[0].avertex().frequencies.b_K1.ws[i];
            IntegranddGammaC_exact_K1<state_datatype> IntegrandK1(Lambda, w, false, Pi);
            state_datatype val_K1 =
                    1. / (2 * M_PI) * integrator<state_datatype, 0>(IntegrandK1, -vmax, vmax, std::std::abs(w / 2), {}, Delta);
            dGammaC_exact.vertex[0].avertex().K1_setvert(0, i, 0, val_K1);
        }


#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta)
        for (int iflat = 0; iflat < (nBOS2 - 2) * (nFER2 - 2); iflat++) {
            int i = 1 + iflat / (nFER2 - 2);
            int j = 1 + iflat - (i - 1) * (nFER2 - 2);
            //for (int i = 1; i<nBOS2-1; i++) {
            //    for (int j = 1; j<nFER2-1; j++) {
            double w = dGammaC_exact.vertex[0].avertex().frequencies.b_K2.ws[i];
            double v = dGammaC_exact.vertex[0].avertex().frequencies.f_K2.ws[j];
            IntegranddGammaC_exact_K2<state_datatype> IntegrandK2(Lambda, w, v, false, Pi);
            state_datatype val_K2 =
                    1. / (2 * M_PI) * integrator<state_datatype, 1>(IntegrandK2, -vmax, vmax, std::std::abs(w / 2), {v}, Delta);
            dGammaC_exact.vertex[0].avertex().K2_setvert(0, i, j, 0, val_K2);
            //    }
        }

#if MAX_DIAG_CLASS>2
#pragma omp parallel for schedule(dynamic) default(none) shared(dGammaC_exact, vmax, Delta)
        for (int iflat = 0; iflat < (nBOS3 - 2) * (nFER3 - 2) * (nFER3 - 2); iflat++) {
            int i = 1 + iflat / (nFER3 - 2) / (nFER3 - 2);
            int j = 1 + iflat / (nFER3 - 2) - (i - 1) * (nFER3 - 2);
            int k = 1 + iflat - (i - 1) * (nFER3 - 2) * (nFER3 - 2) - (j - 1) * (nFER3 - 2);
            //for (int i = 1; i<nBOS3-1; i++) {
            //    for (int j = 1; j<nFER3-1; j++) {
            //        for (int k = 1; k<nFER3-1; k++) {
            double w = dGammaC_exact.vertex[0].avertex().frequencies.b_K3.ws[i];
            double v = dGammaC_exact.vertex[0].avertex().frequencies.f_K3.ws[j];
            double vp = dGammaC_exact.vertex[0].avertex().frequencies.f_K3.ws[k];
            IntegranddGammaC_exact_K3<state_datatype> IntegrandK3(Lambda, w, v, vp, false, Pi);
            state_datatype val_K3 = 1. / (2 * M_PI) *
                                    integrator<state_datatype, 3>(IntegrandK3, -vmax, vmax, std::std::abs(w / 2),
                                                                  {v, vp, std::std::abs(v) - std::std::abs(vp)}, Delta);
            dGammaC_exact.vertex[0].avertex().K3_setvert(0, i, j, k, 0, val_K3);
            //    }
            //}
        }
#endif
        write_hdf("dGammaC_l_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_exact", Lambda, 1,
                  dGammaC_exact);

        State<state_datatype> dGammaC_diff = dGammaC_l - dGammaC_exact;        // final result: contains K1, K2 and K3
        write_hdf("dGammaC_l_version1_U" + std::to_string(glb_U / ((glb_Gamma + Lambda) / 2.)) + ".h5_diff", Lambda, 1,
                  dGammaC_diff);
    }
}

#endif

#endif

#ifdef STATIC_FEEDBACK
/**
 * Compute the right hand side of the flow equations according to Severin Jakobs' channel decomposition with
 * approximated channel feedback and modified self-energy feedback (only static level shift to avoid
 * overbroadening of spectral features)
 * @param Psi    : state at which to compute right hand side
 * @param Lambda : Lambda at which to compute right hand side
 * @return       : dPsi (right hand side of flow equation)
 */
auto rhs_channel_decomposition(const State<Q>& Psi, const double Lambda) -> State<Q> {
    State<Q> dPsi; // result

    SelfEnergy<Q> selfEnergy;
    comp static_shift = real(Psi.selfenergy.valsmooth(0, glb_mu, 0));  // only use a static level shift as self-energy
    selfEnergy.initialize(static_shift, 0.);

    Propagator<Q> G(Lambda, selfEnergy, 'g');    //Initialization of Propagator objects
    Propagator<Q> S(Lambda, selfEnergy, 's');    //Initialization of Propagator objects

    // Self-energy flow
    loop(dPsi.selfenergy, Psi.vertex, S, true);  // self-energy loop

    // Vertex flow
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'a', true); // diff. bubble in the a-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 'p', true); // diff. bubble in the p-channel
    bubble_function(dPsi.vertex, Psi.vertex, Psi.vertex, G, S, 't', true); // diff. bubble in the t-channel

    return dPsi;
}

/**
 * FRG flow according to Severin Jakobs' channel decomposition with approximated channel feedback
 * and modified self-energy feedback (only static level shift to avoid overbroadening of spectral features).
 * Only correct if parameter STATIC_FEEDBACK is defined.
 * @param N_ODE : number of Runge-Kutta ODE iterations
 */
void test_channel_decomposition(int N_ODE) {
    State<Q> state_ini, state_fin;   // create initial and final state
    state_ini.initialize();             // initialize initial state

//    sopt_state(state_ini, Lambda_ini);  // set initial state to SOPT (necessary if Lambda_ini is too small)

    ODE_solver_RK4(state_fin, Lambda_fin, state_ini, Lambda_ini, rhs_channel_decomposition,
                    log_substitution, log_resubstitution, N_ODE); // compute flow

    std::string name = "test_channel_decomposition.h5";
    write_h5_rvecs(name, {"v", "Sigma_re", "Sigma_im", "Sigma_ini_re", "Sigma_ini_im"},
                         {ffreqs,
                          state_fin.selfenergy.Sigma.real(), state_fin.selfenergy.Sigma.imag(),
                          state_ini.selfenergy.Sigma.real(), state_ini.selfenergy.Sigma.imag()});

    write_hdf("channel_decomposition.h5", 0, 1, state_fin);
}
#endif

#endif //KELDYSH_MFRG_TESTFUNCTIONS_H
